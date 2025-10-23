#include "evolve.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdexcept>

#include "io.hpp"

#include "../core/io_utils.hpp"
#include "../core/math_utils.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);
}

// Runge–Kutta 4th order time propagation (ψ' = -i H ψ)
void evolve_rk4_tridiag(const Tridiag& T,
                        const SpectralData& spectral,
                        const Eigen::VectorXcd& psi_init,
                        double dx,
                        double dt,
                        int nsteps,
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
    Eigen::VectorXcd k1(psi.size());
    Eigen::VectorXcd k2(psi.size());
    Eigen::VectorXcd k3(psi.size());
    Eigen::VectorXcd k4(psi.size());
    Eigen::VectorXcd tmp(psi.size());

    double t = 0.0;

    const double norm_sq0 = l2_norm_sq(psi_init, dx);
    const double norm0 = std::sqrt(std::max(0.0, norm_sq0));

    for (int step = 0; step <= nsteps; ++step) {
        if (log_every > 0 && step % log_every == 0) {
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

        tridiag_mul(T, psi, k1);
        k1 *= -I;

        tmp = psi + (0.5 * dt) * k1;
        tridiag_mul(T, tmp, k2);
        k2 *= -I;

        tmp = psi + (0.5 * dt) * k2;
        tridiag_mul(T, tmp, k3);
        k3 *= -I;

        tmp = psi + dt * k3;
        tridiag_mul(T, tmp, k4);
        k4 *= -I;

        Eigen::VectorXcd sum = k1;
        sum.noalias() += 2.0 * k2;
        sum.noalias() += 2.0 * k3;
        sum.noalias() += k4;
        psi.noalias() += (dt / 6.0) * sum;

        double norm = l2_norm(psi, dx);
        if (norm > 0.0) {
            psi /= norm;
        }

        t += dt;
    }
}
