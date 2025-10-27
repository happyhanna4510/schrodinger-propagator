#include "runtime_evolution.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "core/math_utils.hpp"
#include "core/spectral.hpp"
#include "evolve/evolve_factory.hpp"
#include "hamiltonian.hpp"
#include "initial.hpp"
#include "io.hpp"
#include "test_time_reversibility.hpp"

namespace fs = std::filesystem;

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

    Eigen::VectorXcd psi_init = gaussian_on_inner(g, 1.0, 0.35).cast<std::complex<double>>();

    int    K      = (P.K > 0 ? P.K : 4);
    double dt     = (P.dt > 0 ? P.dt : 1e-6);
    double tmax   = (P.tmax > 0 ? P.tmax : 1e-3);
    double tol    = (P.tol > 0 ? P.tol : 1e-12);
    int    nsteps = static_cast<int>(std::llround(tmax / dt));
    const std::string method = P.evolve_method.empty() ? std::string("taylor") : P.evolve_method;

    fs::create_directories(out_dir);
    fs::path csv_path = P.csv_name.empty()
        ? io::make_csv_path(out_dir, P, method, K, dt, g)
        : out_dir / P.csv_name;

    std::vector<double> x_inner;
    const std::vector<double>* x_ptr = nullptr;
    if (P.wide) {
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

#if 0
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
#endif

    if (!P.quiet) {
        std::string method_lower = method;
        std::transform(method_lower.begin(), method_lower.end(), method_lower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        const bool is_taylor = (method_lower == "taylor");
        const bool is_cheb = (method_lower == "cheb" || method_lower == "chebyshev");

        if (is_taylor) {
            std::cout << "\n# === Time evolution (Taylor K=" << K
                      << ", dt=" << dt << ", tmax=" << tmax << ") ===\n";
        } else if (is_cheb) {
            std::cout << "\n# === Time evolution (cheb, dt=" << dt
                      << ", tmax=" << tmax << ", tol=" << tol<<") ===\n";
        } else {
            std::cout << "\n# === Time evolution (" << method
                      << ", dt=" << dt << ", tmax=" << tmax << ") ===\n";
        }
    }

    evolve(method, T, spectral, psi_init, g.dx, dt, tol, nsteps, K,
           csv_path.string(), x_ptr,
           P.wide_re, P.wide_im, P.quiet,
           P.log_every, P.csv_every, P.aggregate,
           P.flush_every, P.no_theta);

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

