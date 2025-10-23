#include "evolve/evolve_factory.hpp"

#include <stdexcept>

#include "evolve/evolve_rk4.hpp"
#include "evolve/evolve_taylor.hpp"

void evolve(const std::string& method,
            const Tridiag& T,
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
    if (method == "taylor") {
        evolve_taylor_tridiag(T, spectral, psi_init, dx, dt, nsteps, K,
                               log_every, csv_path, x_inner,
                               wide_re, wide_im, extras, quiet);
        return;
    }

    if (method == "rk4") {
        evolve_rk4_tridiag(T, spectral, psi_init, dx, dt, nsteps,
                           log_every, csv_path, x_inner,
                           wide_re, wide_im, extras, quiet);
        return;
    }

    if (method == "chebyshev") {
        throw std::runtime_error("chebyshev evolution is not implemented yet");
    }

    throw std::runtime_error("unknown evolution method: " + method);
}

