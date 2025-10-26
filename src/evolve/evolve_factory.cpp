#include "evolve/evolve_factory.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <fstream>
#include <memory>
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
    agg.sum_matvecs += matvecs;
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

bool should_tick(int every, int step, int total_steps) {
    if (every <= 0) {
        return false;
    }
    if (step == 1 || step == total_steps) {
        return true;
    }
    return (step % every) == 0;
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
            bool aggregate) {
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

    std::ofstream csv;
    if (!csv_path.empty()) {
        csv.open(csv_path, std::ios::out | std::ios::trunc);
        if (!csv) {
            throw std::runtime_error("failed to open log csv");
        }
        write_step_csv_header(csv, is_cheb);
    }

    io::WideDump wide(csv_path, x_inner, wide_re, wide_im);

    Eigen::VectorXcd psi = psi_init;
    double t = 0.0;

    IntervalAgg console_agg;
    IntervalAgg csv_agg;

    for (int step = 1; step <= nsteps; ++step) {
        const auto start = std::chrono::steady_clock::now();
        StepResult result = evolver->step(psi);
        const auto end = std::chrono::steady_clock::now();

        t += dt;
        const double dt_ms = std::chrono::duration<double, std::milli>(end - start).count();

        const double norm_sq = l2_norm_sq(psi, dx);
        const double norm_err = std::abs(norm_sq - 1.0);

        update_interval(console_agg, dt_ms, result.matvecs, norm_err, result);
        if (csv.is_open()) {
            update_interval(csv_agg, dt_ms, result.matvecs, norm_err, result);
        }

        const bool console_tick = should_tick(cfg.log_every, step, nsteps);
        const bool csv_tick = csv.is_open() && should_tick(cfg.csv_every, step, nsteps);
        const bool need_metrics = csv_tick || (!quiet && console_tick);

        StepMetrics metrics_current{};
        if (need_metrics) {
            metrics_current = compute_step_metrics(spectral, psi, t, dx);
        }

        if (csv_tick) {
            IntervalStats stats = finalize_interval(csv_agg, cfg.aggregate);
            StepMetrics metrics_csv = metrics_current;
            metrics_csv.norm_err = stats.norm_err;
            write_step_csv_row(csv, method_norm, step, t, dt, stats.dt_ms,
                               stats.matvecs, metrics_csv, stats.K_used, stats.bn_ratio, is_cheb);
            reset_interval(csv_agg);
        }

        if (console_tick) {
            IntervalStats stats = finalize_interval(console_agg, cfg.aggregate);
            if (!quiet) {
                StepMetrics metrics_console = metrics_current;
                metrics_console.norm_err = stats.norm_err;
                print_step_console(method_norm, step, t, stats.dt_ms, stats.matvecs,
                                   metrics_console, stats.K_used);
            }
            reset_interval(console_agg);
        }

        if (wide.enabled()) {
            wide.write(psi, t);
        }
    }
}

