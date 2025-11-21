#include "evolve/evolve_factory.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <filesystem>

#include "core/io_utils.hpp"
#include "core/math_utils.hpp"
#include "evolve/evolve_chebyshev.hpp"
#include "evolve/evolve_rk4.hpp"
#include "evolve/evolve_taylor.hpp"
#include "io.hpp"

namespace fs = std::filesystem;


namespace {

struct IntervalAgg {
    int count = 0;
    double sum_ms = 0.0;
    long long sum_matvecs = 0;
    double max_norm_err = 0.0;
    double last_dt_ms = 0.0;
    double last_matvecs = 0.0;
    double last_norm_err = 0.0;
    std::optional<int> last_K_used;
    std::optional<double> last_bn_ratio;
};

struct IntervalStats {
    double dt_ms = 0.0;
    double matvecs = 0.0;
    double norm_err = 0.0;
    std::optional<int> K_used;
    std::optional<double> bn_ratio;
};

void reset_interval(IntervalAgg& agg) {
    agg = IntervalAgg{};
}

void update_interval(IntervalAgg& agg,
                     double dt_ms,
                     int matvecs,
                     double norm_err,
                     const StepResult& result) {
    ++agg.count;
    agg.sum_ms += dt_ms;
    agg.sum_matvecs += static_cast<long long>(matvecs);
    agg.max_norm_err = std::max(agg.max_norm_err, norm_err);
    agg.last_dt_ms = dt_ms;
    agg.last_matvecs = static_cast<double>(matvecs);
    agg.last_norm_err = norm_err;
    agg.last_K_used = result.K_used;
    agg.last_bn_ratio = result.bn_ratio;
}

IntervalStats finalize_interval(const IntervalAgg& agg, bool aggregate) {
    IntervalStats stats;
    if (agg.count > 0) {
        if (aggregate) {
            stats.dt_ms = agg.sum_ms / static_cast<double>(agg.count);
            stats.matvecs = static_cast<double>(agg.sum_matvecs) / static_cast<double>(agg.count);
            stats.norm_err = agg.max_norm_err;
        } else {
            stats.dt_ms = agg.last_dt_ms;
            stats.matvecs = agg.last_matvecs;
            stats.norm_err = agg.last_norm_err;
        }
    }
    stats.K_used = agg.last_K_used;
    stats.bn_ratio = agg.last_bn_ratio;
    return stats;
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

std::string normalize_method(const std::string& method) {
    const std::string lower = to_lower(method);
    if (lower == "chebyshev") {
        return "cheb";
    }
    return lower;
}

} // namespace

void evolve(const std::string& method,
            const Tridiag& T,
            const SpectralData& spectral,
            const Eigen::VectorXcd& psi_init,
            double dx,
            double dt,
            double tol,
            int nsteps,
            int K,
            const std::string& csv_path,
            const std::vector<double>* x_inner,
            bool wide_re,
            bool wide_im,
            bool quiet,
            int log_every,
            int csv_every,
            bool aggregate,
            int flush_every,
            bool no_theta,
            bool profile,
            const EnergyLogConfig& energy_cfg,
            const DensityLogConfig& density_cfg) {
    const std::string method_norm = normalize_method(method);
    const bool is_cheb = (method_norm == "cheb");

    EvolverConfig cfg;
    cfg.dt = dt;
    cfg.dx = dx;
    cfg.K = K;
    cfg.hbar = 1.0;
    cfg.tolerance = tol;
    cfg.log_every = log_every;
    cfg.csv_every = csv_every;
    cfg.aggregate = aggregate;
    cfg.flush_every = flush_every;
    cfg.no_theta = no_theta;
    cfg.profile = profile;

    const char* theta_debug_env = std::getenv("THETA_DEBUG");
    const bool theta_debug_enabled = theta_debug_env && theta_debug_env[0] != '\0';
    double theta_debug_lo = 3.5;
    double theta_debug_hi = 4.5;
    if (const char* win_env = std::getenv("THETA_DEBUG_WINDOW")) {
        std::string win_str(win_env);
        const auto pos = win_str.find_first_of(",;:");
        if (pos != std::string::npos) {
            try {
                theta_debug_lo = std::stod(win_str.substr(0, pos));
                theta_debug_hi = std::stod(win_str.substr(pos + 1));
                if (theta_debug_lo > theta_debug_hi) {
                    std::swap(theta_debug_lo, theta_debug_hi);
                }
            } catch (...) {
                std::cerr << "# [theta-debug] failed to parse THETA_DEBUG_WINDOW='"
                          << win_str << "' (expected <lo>,<hi>)\n";
            }
        } else {
            std::cerr << "# [theta-debug] THETA_DEBUG_WINDOW='" << win_str
                      << "' missing delimiter ',' ';' or ':'\n";
        }
    }

    std::unique_ptr<EvolverBase> evolver;

    if (method_norm == "taylor") {
        evolver = std::make_unique<TaylorEvolver>(T, cfg);
    } else if (method_norm == "rk4") {
        evolver = std::make_unique<Rk4Evolver>(T, cfg);
    } else if (method_norm == "cheb") {
        evolver = std::make_unique<ChebyshevEvolver>(T, spectral, cfg);
    } else {
        throw std::runtime_error("unknown evolution method: " + method);
    }

    const bool enable_csv = (!csv_path.empty() && csv_every != 0);

    std::ofstream csv;
    int csv_rows = 0;
    if (enable_csv) {
        csv.open(csv_path, std::ios::out | std::ios::trunc);
        if (!csv) {
            throw std::runtime_error("failed to open log csv");
        }
        write_step_csv_header(csv, is_cheb);
    }

    io::WideDump wide(csv_path, x_inner, wide_re, wide_im);

    const bool density_logging = density_cfg.enabled &&
                                 !density_cfg.x_inner.empty() &&
                                 !density_cfg.num_csv_path.empty() &&
                                 !density_cfg.ref_csv_path.empty();
    const int density_every = density_logging ? std::max(1, density_cfg.every) : 1;

    std::ofstream density_num_csv;
    std::ofstream density_ref_csv;
    bool density_header_written = false;
    int density_rows = 0;

    if (density_logging) {
        if (static_cast<std::size_t>(psi_init.size()) != density_cfg.x_inner.size()) {
            throw std::runtime_error("density logging grid mismatch");
        }

        fs::path num_path = fs::path(density_cfg.num_csv_path);
        fs::path ref_path = fs::path(density_cfg.ref_csv_path);
        std::error_code ec;
        if (!num_path.parent_path().empty()) {
            fs::create_directories(num_path.parent_path(), ec);
        }
        if (!ref_path.parent_path().empty()) {
            fs::create_directories(ref_path.parent_path(), ec);
        }

        density_num_csv.open(num_path, std::ios::out | std::ios::trunc);
        density_ref_csv.open(ref_path, std::ios::out | std::ios::trunc);
        if (!density_num_csv || !density_ref_csv) {
            throw std::runtime_error("failed to open density CSV outputs");
        }
        density_num_csv << std::setprecision(16);
        density_ref_csv << std::setprecision(16);
    }

    Eigen::VectorXcd psi = psi_init;
    double t = 0.0;

    const bool energy_logging = energy_cfg.enabled && energy_cfg.potential && !energy_cfg.csv_path.empty();
    std::ofstream energy_csv;
    int energy_rows = 0;
    fs::path energy_path;
    if (energy_logging) {
        energy_path = fs::path(energy_cfg.csv_path);
        std::error_code ec;
        if (!energy_path.parent_path().empty()) {
            fs::create_directories(energy_path.parent_path(), ec);
        }
        energy_csv.open(energy_path, std::ios::out | std::ios::trunc);
        if (!energy_csv) {
            throw std::runtime_error("failed to open energy log csv");
        }
        energy_csv << "t,E\n";
        energy_csv << std::setprecision(16);
    }

    if (energy_logging && energy_cfg.potential->size() != static_cast<std::size_t>(psi.size())) {
        throw std::runtime_error("energy logging potential size mismatch");
    }

    IntervalAgg agg;
    int tick_counter = 0;

    struct ProfileAccumulator {
        double step_us = 0.0;
        double work_us = 0.0;
        double rhs_us = 0.0;
        long long reallocations = 0;
        long long samples = 0;
    } profile_acc;

    struct ThetaDebugSummary {
        double min_theta_abs = std::numeric_limits<double>::infinity();
        double max_dev_z_coeff = 0.0;
        double max_dev_z_grid = 0.0;
        bool anomalies_z_grid = false;
        bool theta_consistent_with_coeff = true;
        long long samples = 0;
    } theta_summary;


    auto compute_reference_state = [&](double time) {
        constexpr std::complex<double> I(0.0, 1.0);
        Eigen::ArrayXcd phase = (-I * spectral.evals.array() * time).exp();
        Eigen::VectorXcd coeffs = (spectral.coeffs0.array() * phase).matrix();
        return (spectral.eigenvectors * coeffs).eval();
    };


    auto write_density_row = [&](const Eigen::VectorXcd& psi, const Eigen::VectorXcd& psi_ref, double time) {
        if (!density_header_written) {
            auto write_header = [&](std::ofstream& f) {
                f << "t";
                for (double xi : density_cfg.x_inner) {
                    f << ',' << xi;
                }
                f << '\n';
            };
            write_header(density_num_csv);
            write_header(density_ref_csv);
            density_header_written = true;
        }

        auto write_row = [&](std::ofstream& f, const Eigen::VectorXcd& state) {
            f << time;
            for (int i = 0; i < state.size(); ++i) {
                const double re = state[i].real();
                const double im = state[i].imag();
                f << ',' << (re * re + im * im);
            }
            f << '\n';
        };

        write_row(density_num_csv, psi);
        write_row(density_ref_csv, psi_ref);
        ++density_rows;
        if (flush_every > 0 && (density_rows % flush_every) == 0) {
            density_num_csv.flush();
            density_ref_csv.flush();
        }
    };

    for (int step = 0; step < nsteps; ++step)
    {
        const auto start = std::chrono::steady_clock::now();
        StepResult result = evolver->step(psi);
        const auto end = std::chrono::steady_clock::now();

        t += dt;
        const double dt_ms = std::chrono::duration<double, std::milli>(end - start).count();

        const double norm_sq = l2_norm_sq(psi, dx);
        const double norm_err = std::abs(norm_sq - 1.0);

        update_interval(agg, dt_ms, result.matvecs, norm_err, result);

        const bool log_density = density_logging &&
                                 ((step % density_cfg.every) == 0 || step == 0 || step + 1 == nsteps);
        if (log_density) {
            Eigen::VectorXcd psi_ref = compute_reference_state(t);
            write_density_row(psi, psi_ref, t);
        }


        if (cfg.profile && result.profile) {
            profile_acc.step_us += result.profile->step_us;
            profile_acc.work_us += result.profile->work_us;
            profile_acc.rhs_us += result.profile->rhs_us;
            profile_acc.reallocations += result.profile->reallocations;
            ++profile_acc.samples;
        }

        const bool tick = (cfg.log_every > 0) &&
                          ((step % cfg.log_every) == 0 || step == 0 || step + 1 == nsteps);


        if (tick)
        {
            IntervalStats stats = finalize_interval(agg, cfg.aggregate);
            double theta_rel = std::numeric_limits<double>::quiet_NaN();
            double theta_abs = std::numeric_limits<double>::quiet_NaN();
            std::optional<double> e_true;
            if (!cfg.no_theta) {
                theta_rel = compute_theta(spectral, psi, t, dx, /*relative=*/true);
                theta_abs = compute_theta(spectral, psi, t, dx, /*relative=*/false);

            }
            const int step_out = step + 1;

            
            if (theta_debug_enabled && std::isfinite(theta_abs)) {
                const double t_num = t;
                const double t_ref = static_cast<double>(step_out) * dt;
                std::streamsize old_prec = std::cout.precision();
                std::ios::fmtflags old_flags = std::cout.flags();
                std::cout.setf(std::ios::scientific);
                std::cout << "# [theta-debug-time] step=" << step_out
                          << " t_num=" << std::setprecision(15) << t_num
                          << " t_ref=" << std::setprecision(15) << t_ref
                          << " theta_abs=" << std::setprecision(6) << theta_abs
                          << '\n';
                std::cout.precision(old_prec);
                std::cout.flags(old_flags);
            }

            const bool in_theta_window = theta_debug_enabled &&
                                         (t >= theta_debug_lo) && (t <= theta_debug_hi);
            if (in_theta_window) {
                Eigen::VectorXcd psi_ref = compute_reference_state(t);
                std::complex<double> z_grid = inner_dx(psi_ref, psi, dx);
                std::complex<double> z_coeff = compute_overlap_coeffs(spectral, psi, t, dx);
                const double theta_raw_grid = std::arg(z_grid);
                const double theta_raw_coeff = std::arg(z_coeff);
                const double z_abs_grid = std::abs(z_grid);
                const double z_abs_coeff = std::abs(z_coeff);
                const double norm_psi = l2_norm(psi, dx);
                const double norm_ref = l2_norm(psi_ref, dx);
                theta_summary.min_theta_abs = std::min(theta_summary.min_theta_abs, theta_abs);
                theta_summary.max_dev_z_coeff = std::max(theta_summary.max_dev_z_coeff, std::abs(z_abs_coeff - 1.0));
                theta_summary.max_dev_z_grid = std::max(theta_summary.max_dev_z_grid, std::abs(z_abs_grid - 1.0));
                theta_summary.anomalies_z_grid = theta_summary.anomalies_z_grid || (z_abs_grid < 0.9 || z_abs_grid > 1.1);
                theta_summary.theta_consistent_with_coeff = theta_summary.theta_consistent_with_coeff &&
                                                           !((z_abs_coeff > 0.9) && (theta_abs > 1e-3));
                ++theta_summary.samples;
                std::streamsize old_prec = std::cout.precision();
                std::ios::fmtflags old_flags = std::cout.flags();
                std::cout.setf(std::ios::scientific);
                std::cout << "# [theta-debug-overlap] t=" << std::setprecision(9) << t
                          << " |z_grid|=" << std::setprecision(6) << z_abs_grid
                          << " |z_coeff|=" << std::setprecision(6) << z_abs_coeff
                          << " theta_raw_grid=" << std::setprecision(6) << theta_raw_grid
                          << " theta_raw_coeff=" << std::setprecision(6) << theta_raw_coeff
                          << " norm_psi=" << std::setprecision(6) << norm_psi
                          << " norm_ref=" << std::setprecision(6) << norm_ref
                          << " theta_abs=" << std::setprecision(6) << theta_abs
                          << '\n';
                std::cout.precision(old_prec);
                std::cout.flags(old_flags);
            }


            if (csv.is_open() && ((cfg.csv_every <= 1) || (tick_counter % cfg.csv_every) == 0)) {
                write_step_csv_row(csv, method_norm, step_out, t, dt, stats.dt_ms,
                                   stats.matvecs, stats.norm_err, theta_rel, theta_abs,
                                   stats.K_used, stats.bn_ratio, is_cheb);
                ++csv_rows;
                if (cfg.flush_every > 0 && (csv_rows % cfg.flush_every) == 0) {
                    csv.flush();
                }
            }

            if (energy_logging && energy_csv.is_open()) {
                const double energy_val = compute_energy(psi, *energy_cfg.potential,
                                                         dx, energy_cfg.hbar, energy_cfg.mass);
                energy_csv << t << ',' << energy_val << '\n';
                ++energy_rows;
                if (cfg.flush_every > 0 && (energy_rows % cfg.flush_every) == 0) {
                    energy_csv.flush();
                }
            }

            if (!quiet && cfg.log_every > 0) {
                print_step_console(method_norm, step_out, t, stats.dt_ms, stats.matvecs,
                                   stats.norm_err, theta_rel, theta_abs, stats.K_used);
            }
            if (wide.enabled()) {
                wide.write(psi, t);
            }

            reset_interval(agg);
            ++tick_counter;
        }


    }

    if (theta_debug_enabled && theta_summary.samples > 0) {
        const bool coeff_stable = theta_summary.max_dev_z_coeff < 1e-2;
        const bool anomalies_grid = theta_summary.anomalies_z_grid;
        const bool theta_consistent = theta_summary.theta_consistent_with_coeff;
        const double min_theta_abs = theta_summary.min_theta_abs;

        std::string probable_cause;
        if (anomalies_grid && coeff_stable) {
            probable_cause = "grid overlap sensitivity to dx normalization causing suppressed |z_grid|";
        } else if (anomalies_grid) {
            probable_cause = "grid overlap instability (possible amplitude scaling or dx sensitivity)";
        } else {
            probable_cause = "overlaps consistent; no suppression observed";
        }

        std::string recommendation;
        if (anomalies_grid) {
            recommendation = "prefer coefficient-space z(t) for diagnostics; revisit grid normalization and dx scaling";
        } else {
            recommendation = "diagnostics stable; continue using both overlaps for cross-checks";
        }

        std::streamsize old_prec = std::cout.precision();
        std::ios::fmtflags old_flags = std::cout.flags();
        std::cout.setf(std::ios::scientific);
        std::cout << "\n=== FINAL SUMMARY ===\n";
        std::cout << "theta_abs minimal: " << std::setprecision(6) << min_theta_abs << "\n";
        std::cout << "max deviation of |z_coeff| from 1: " << std::setprecision(6)
                  << theta_summary.max_dev_z_coeff << "\n";
        std::cout << "observed anomalies in |z_grid|: " << (anomalies_grid ? "yes" : "no") << "\n";
        std::cout << "|z_coeff| stayed close to 1: " << (coeff_stable ? "yes" : "no") << "\n";
        std::cout << "theta consistent with z_coeff: " << (theta_consistent ? "yes" : "no") << "\n";
        std::cout << "probable cause: " << probable_cause << "\n";
        std::cout << "recommendation: " << recommendation << "\n";
        std::cout << "======================\n";
        std::cout.precision(old_prec);
        std::cout.flags(old_flags);
    }

    if (cfg.profile && profile_acc.samples > 0 && !quiet) {
        const double inv = 1.0 / static_cast<double>(profile_acc.samples);
        std::cout << "\n# Profiling averages (" << profile_acc.samples << " steps)\n";
        std::cout << "#   section        avg time [us]\n";
        std::cout << "#   step_total     " << profile_acc.step_us * inv << "\n";
        std::cout << "#   work_total     " << profile_acc.work_us * inv << "\n";
        std::cout << "#   rhs_matvec     " << profile_acc.rhs_us * inv << "\n";
        std::cout << "#   reallocations  "
                  << (profile_acc.reallocations * inv) << "\n";
    }
}

