#include "evolve.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>

#include "io.hpp"

#include "../core/io_utils.hpp"
#include "../core/math_utils.hpp"
#include "../core/workspace.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);

void taylor_step_tridiag(const Tridiag& T,
                         Eigen::VectorXcd& psi,
                         double dt,
                         int K,
                         TaylorWorkspace& workspace) {
    workspace.resize(psi.size());
    auto& sum = workspace.sum();
    auto& vk  = workspace.vk();
    auto& tmp = workspace.tmp();

    sum = psi;
    vk  = psi;

    const std::complex<double> scale = -I * dt;
    for (int k = 1; k <= K; ++k) {
        tridiag_mul(T, vk, tmp);
        tmp *= (scale / static_cast<double>(k));
        vk = tmp;
        sum.noalias() += vk;
    }

    psi.swap(sum);
}
} // namespace

void evolve_taylor_tridiag(const Tridiag& T,
                           const SpectralData& spectral,
                           const Eigen::VectorXcd& psi_init,
                           double dx,
                           double dt,
                           int nsteps,
                           int K,
                           int log_every,
                           const std::string& csv_path,
                           const std::vector<double>* x_inner,
                           bool wide_re,
                           bool wide_im,
                           LogExtras extras,
                           bool quiet) {
    std::ofstream csv;
    if (!csv_path.empty()) {
        csv.open(csv_path, std::ios::out | std::ios::trunc);
        if (!csv) {
            throw std::runtime_error("failed to open log csv");
        }
        write_csv_header(csv, extras);
    }

    io::WideDump wide(csv_path, x_inner, wide_re, wide_im);

    Eigen::VectorXcd psi = psi_init;
    double t = 0.0;

    TaylorWorkspace workspace;
    workspace.resize(psi.size());

    const double norm_sq0 = l2_norm_sq(psi_init, dx);
    const double norm0 = std::sqrt(std::max(0.0, norm_sq0));

    for (int step = 0; step <= nsteps; ++step) {
        if (step % log_every == 0) {
            LogSnapshot snap = collect_snapshot(spectral, psi, t, dx, step, norm_sq0, norm0);

            if (csv.is_open()) {
                write_csv_row(csv, snap, extras);
            }
            if (!quiet) {
                print_snapshot(snap, extras);
            }
            if (wide.enabled()) {
                wide.write(psi, t);
            }
        }

        if (step == nsteps) {
            break;
        }

        taylor_step_tridiag(T, psi, dt, K, workspace);
        t += dt;
    }
}

RevTestResult test_time_reversibility(const Tridiag& T,
                                      const Eigen::VectorXcd& psi_init,
                                      double dx,
                                      double dt,
                                      int nsteps,
                                      int K) {
    RevTestResult r{};
    r.norm0 = l2_norm(psi_init, dx);

    Eigen::VectorXcd psi = psi_init;
    TaylorWorkspace workspace;
    workspace.resize(psi.size());

    auto update_drift = [&](double& cur_max, const Eigen::VectorXcd& v) {
        double n = l2_norm(v, dx);
        cur_max = std::max(cur_max, std::abs(n - 1.0));
    };

    double max_drift = 0.0;
    for (int i = 0; i < nsteps; ++i) {
        taylor_step_tridiag(T, psi, dt, K, workspace);
        if ((i & 0x3FF) == 0) {
            update_drift(max_drift, psi);
        }
    }
    r.norm_fwd = l2_norm(psi, dx);
    update_drift(max_drift, psi);

    for (int i = 0; i < nsteps; ++i) {
        taylor_step_tridiag(T, psi, -dt, K, workspace);
        if ((i & 0x3FF) == 0) {
            update_drift(max_drift, psi);
        }
    }
    r.norm_back = l2_norm(psi, dx);
    update_drift(max_drift, psi);

    r.err_back = (psi - psi_init).norm();
    r.max_drift = max_drift;
    return r;
}
