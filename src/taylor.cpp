#include "taylor.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "evolve.hpp"
#include "hamiltonian.hpp"
#include "initial.hpp"
#include "io.hpp"

#include "core/math_utils.hpp"
#include "core/spectral.hpp"

namespace fs = std::filesystem;

Tridiag make_tridiag_from_dense(const Eigen::Ref<const Eigen::MatrixXd>& H) {
    const int M = static_cast<int>(H.rows());
    const int off_size = std::max(0, M - 1);
    Tridiag T{Eigen::VectorXd(M), Eigen::VectorXd(off_size), Eigen::VectorXd(off_size)};
    T.a = H.diagonal();
    if (off_size > 0) {
        T.b = H.diagonal(1);
        T.c = H.diagonal(-1);
    }
    return T;
}

fs::path run_time_evolution(const Grid& g,
                            const std::vector<double>& U_true,
                            const Params& P,
                            const fs::path& out_dir) {
    std::vector<double> U_evol = U_true;
    double cur_max = *std::max_element(U_evol.begin(), U_evol.end());
    double scale = (cur_max > 0.0) ? (P.Umax / cur_max) : 1.0;
    for (double& v : U_evol) {
        v *= scale;
    }

    if (!P.quiet) {
        std::cout << "# Evolution scaling: current max(U)=" << cur_max
                  << ", target Umax=" << P.Umax
                  << ", scale=" << scale << "\n";
    }

    Eigen::MatrixXd H_evol = build_hamiltonian(g, U_evol);
    Tridiag T = make_tridiag_from_dense(H_evol);

    Eigen::VectorXcd psi_init = gaussian_on_inner(g, 1.0, 0.35).cast<cplx>();

    int    K      = (P.K > 0 ? P.K : 4);
    double dt     = (P.dt > 0 ? P.dt : 1e-6);
    double tmax   = (P.tmax > 0 ? P.tmax : 1e-3);
    int    nsteps = static_cast<int>(std::llround(tmax / dt));
    int    log_every = (P.log_every > 0 ? P.log_every : 100);
    const std::string method = P.evolve_method.empty() ? std::string("taylor") : P.evolve_method;

    fs::create_directories(out_dir);
    fs::path csv_path = P.csv_name.empty()
        ? io::make_csv_path(out_dir, P, method, K, dt, g)
        : out_dir / P.csv_name;

    std::vector<double> x_inner;
    const std::vector<double>* x_ptr = nullptr;
    if (!P.no_wide) {
        x_inner.resize(g.N - 2);
        for (int i = 0; i < g.N - 2; ++i) {
            x_inner[i] = g.x[i + 1];
        }
        x_ptr = &x_inner;
    }

    SpectralData spectral = make_spectral_data(H_evol, psi_init, g.dx);
    double Emin = spectral.evals.minCoeff();
    double Emax = spectral.evals.maxCoeff();
    double dE   = Emax - Emin;

    if (!P.quiet) {
        std::cout << "# Spectrum: Emin=" << Emin << "  Emax=" << Emax
                  << "  dE=" << dE << "\n";
    }
    maybe_warn_timestep(dt, dE, P.quiet);

    if (method == "taylor") {
        RevTestResult rt = test_time_reversibility(T, psi_init, g.dx, dt, nsteps, K);
        if (!P.quiet) {
            std::cout << "\n# Reversibility test (K=" << K << ", dt=" << dt
                      << ", nsteps=" << nsteps << ")\n"
                      << "  err_back = " << std::scientific << rt.err_back << "\n"
                      << std::defaultfloat
                      << "  norm0    = " << std::setprecision(12) << rt.norm0 << "\n"
                      << "  norm_fwd = " << std::setprecision(12) << rt.norm_fwd << "\n"
                      << "  norm_back= " << std::setprecision(12) << rt.norm_back << "\n"
                      << std::scientific
                      << "  max|norm-1| (both passes) = " << rt.max_drift << "\n";
            std::cout << std::defaultfloat;
        }
    } else if (!P.quiet) {
        std::cout << "\n# Reversibility test is not implemented for method '"
                  << method << "'.\n";
    }

    if (!P.quiet) {
        if (method == "taylor") {
            std::cout << "\n# === Time evolution (Taylor K=" << K
                      << ", dt=" << dt << ", tmax=" << tmax << ") ===\n";
        } else {
            std::cout << "\n# === Time evolution (" << method
                      << ", dt=" << dt << ", tmax=" << tmax << ") ===\n";
        }
    }

    LogExtras extras = LogExtras::None;
    if (P.log_p0) {
        extras = static_cast<LogExtras>(static_cast<unsigned>(extras) |
                                        static_cast<unsigned>(LogExtras::P0));
    }
    if (P.log_err_exact) {
        extras = static_cast<LogExtras>(static_cast<unsigned>(extras) |
                                        static_cast<unsigned>(LogExtras::ErrPhi));
    }

    evolve(method, T, spectral, psi_init, g.dx, dt, nsteps, K,
           log_every, csv_path.string(), x_ptr,
           P.wide_re, P.wide_im, extras, P.quiet);

    if (!P.quiet) {
        std::cout << "# log saved to: " << csv_path.string() << "\n";
    }

    return csv_path;
}

fs::path run_taylor_evolution(const Grid& g,
                              const std::vector<double>& U_true,
                              const Params& P,
                              const fs::path& out_dir) {
    Params P_taylor = P;
    P_taylor.evolve_method = "taylor";
    return run_time_evolution(g, U_true, P_taylor, out_dir);
}
