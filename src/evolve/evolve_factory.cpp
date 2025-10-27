#include "evolve/evolve_factory.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "core/io_utils.hpp"
#include "core/math_utils.hpp"
#include "evolve/evolve_chebyshev.hpp"
#include "evolve/evolve_rk4.hpp"
#include "evolve/evolve_taylor.hpp"
#include "io.hpp"

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
            bool no_theta) {
    const std::string method_norm = normalize_method(method);
    const bool is_cheb = (method_norm == "cheb");

    EvolverConfig cfg;
    cfg.dt = dt;
    cfg.dx = dx;
    cfg.K = K;
    cfg.hbar = 1.0;
    cfg.tolerance = 1e-12;
    cfg.log_every = log_every;
    cfg.csv_every = csv_every;
    cfg.aggregate = aggregate;
    cfg.flush_every = flush_every;
    cfg.no_theta = no_theta;

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

    Eigen::VectorXcd psi = psi_init;
    double t = 0.0;

    IntervalAgg agg;
    int tick_counter = 0;

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

        const bool tick = (cfg.log_every > 0) &&
                          ((step % cfg.log_every) == 0 || step == 0 || step + 1 == nsteps);


        if (tick) 
        {
            IntervalStats stats = finalize_interval(agg, cfg.aggregate);
            double theta = std::numeric_limits<double>::quiet_NaN();
            if (!cfg.no_theta) {
                theta = compute_theta(spectral, psi, t, dx);
            }
            const int step_out = step + 1;

            if (csv.is_open() && ((cfg.csv_every <= 1) || (tick_counter % cfg.csv_every) == 0)) {
                write_step_csv_row(csv, method_norm, step_out, t, dt, stats.dt_ms,
                                   stats.matvecs, stats.norm_err, theta,
                                   stats.K_used, stats.bn_ratio, is_cheb);
                ++csv_rows;
                if (cfg.flush_every > 0 && (csv_rows % cfg.flush_every) == 0) {
                    csv.flush();
                }
            }

            if (!quiet && cfg.log_every > 0) {
                print_step_console(method_norm, step_out, t, stats.dt_ms, stats.matvecs,
                                   stats.norm_err, theta, stats.K_used);
            }
            if (wide.enabled()) {
                wide.write(psi, t);
            }

            reset_interval(agg);
            ++tick_counter;
        }


    }
}

