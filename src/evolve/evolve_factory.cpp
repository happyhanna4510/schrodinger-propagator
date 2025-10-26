#include "evolve/evolve_factory.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

#include "core/io_utils.hpp"
#include "evolve/evolve_chebyshev.hpp"
#include "evolve/evolve_rk4.hpp"
#include "evolve/evolve_taylor.hpp"
#include "io.hpp"

namespace {

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
            bool quiet) {
    const std::string method_norm = normalize_method(method);
    const bool is_cheb = (method_norm == "cheb");

    EvolverConfig cfg;
    cfg.dt = dt;
    cfg.dx = dx;
    cfg.K = K;
    cfg.hbar = 1.0;
    cfg.tolerance = 1e-12;

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

    for (int step = 1; step <= nsteps; ++step) {
        const auto start = std::chrono::steady_clock::now();
        StepResult result = evolver->step(psi);
        const auto end = std::chrono::steady_clock::now();

        t += dt;
        const double dt_ms = std::chrono::duration<double, std::milli>(end - start).count();

        StepMetrics metrics = compute_step_metrics(spectral, psi, t, dx);

        if (csv.is_open()) {
            write_step_csv_row(csv, method_norm, step, t, dt, dt_ms, result.matvecs,
                               metrics, result.K_used, result.bn_ratio, is_cheb);
        }

        if (!quiet) {
            print_step_console(method_norm, step, t, dt_ms, result.matvecs,
                               metrics, result.K_used);
        }

        if (wide.enabled()) {
            wide.write(psi, t);
        }
    }
}

