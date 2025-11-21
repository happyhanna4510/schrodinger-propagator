#include "runtime_evolution.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues> // [ TEST] required for eigenstate-based initialization

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept> // [ TEST] for eigen init error handling
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
    const auto minmax_true = std::minmax_element(U_true.begin(), U_true.end());
    double cur_max = 0.0;
    if (!U_true.empty()) {
        cur_max = *minmax_true.second;
    }

    double scale = 1.0;
    bool scaling_applied = false;
    if (P.Umax_specified && P.Umax > 0.0) {
        // [AI PATCH] Only scale the Morse potential when the user requested a cap via --Umax/--Vcap.
        scale = (cur_max > 0.0) ? (P.Umax / cur_max) : 1.0;
        for (double& v : U_evol) {
            v *= scale;
        }
        scaling_applied = true;
    }

    const auto minmax_scaled = std::minmax_element(U_evol.begin(), U_evol.end());
    double scaled_min = 0.0;
    double scaled_max = 0.0;
    if (!U_evol.empty()) {
        scaled_min = *minmax_scaled.first;
        scaled_max = *minmax_scaled.second;
    }

    if (!P.quiet) {
        if (scaling_applied) {
            std::cout << "# [PATCH] Evolution scaling: current max(U)=" << cur_max
                      << ", target Umax=" << P.Umax
                      << ", scale=" << scale << "\n";
        } else {
            std::cout << "# [PATCH] Evolution scaling: using raw Morse potential (scale=1)\n";
        }
        // [ PATCH] added diagnostic to confirm the Morse potential is active during evolution
        std::cout << "# [PATCH] [evolve] using Morse potential with gamma=" << P.gamma
                  << ", V_min=" << scaled_min
                  << ", V_max=" << scaled_max << "\n";
    }

    const double field_coeff = (g.xmax != 0.0) ? (P.U0 / g.xmax) : 0.0;
    std::vector<double> U_energy_inner;
    U_energy_inner.reserve(std::max(0, g.N - 2));
    for (int i = 1; i < g.N - 1; ++i) {
        const std::size_t idx = static_cast<std::size_t>(i);
        U_energy_inner.push_back(U_evol[idx] + field_coeff * g.x[idx]);
    }

    Eigen::MatrixXd H_evol = build_hamiltonian(g, U_evol, P.U0);
    Tridiag T = make_tridiag_from_dense(H_evol);

    Eigen::VectorXcd psi_init;
    if (P.init == "real-gauss") {
        psi_init = gaussian_on_inner(g, P.x0, P.sigma).cast<std::complex<double>>();
    } else if (P.init == "eigen") {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_init(H_evol); // [ TEST] diagonalize for eigen init
        if (es_init.info() != Eigen::Success) {
            throw std::runtime_error("eigensolve failed while preparing eigen init"); // [TEST] propagate solver issues
        }
        Eigen::MatrixXd Vreal = es_init.eigenvectors();
        renormalize_eigenvectors(Vreal, g.dx); // [ TEST] ensure orthonormality on the grid
        int n = P.eigen_index;
        if (n < 0 || n >= Vreal.cols()) {
            throw std::runtime_error("--eigen_index out of range for available Morse eigenstates"); // [TEST] guard invalid index
        }
        psi_init = Vreal.col(n).cast<std::complex<double>>();
        double norm = psi_init.norm();
        if (norm > 0.0) {
            psi_init /= norm; // [ TEST] safety normalization
        }
        if (!P.quiet) {
            double energy = es_init.eigenvalues()(n);
            std::cout << "# [PATCH TEST] [init] using Morse eigenstate n=" << n
                      << ", E_num=" << energy << "\n";
        }
    } else {
        psi_init = gaussian_complex_on_inner(g, P.x0, P.sigma, P.k0);
    }

    
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
    if (P.wide || P.export_ref_density) {
        x_inner.resize(g.N - 2);
        for (int i = 0; i < g.N - 2; ++i) {
            x_inner[i] = g.x[i + 1];
        }
        if (P.wide) {
            x_ptr = &x_inner;
        }
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

    EnergyLogConfig energy_cfg;
    if (P.log_energy) {
        energy_cfg.enabled = true;
        energy_cfg.potential = &U_energy_inner;
        energy_cfg.hbar = 1.0;
        energy_cfg.mass = 1.0;
        const fs::path energy_path = io::make_energy_csv_path(csv_path);
        energy_cfg.csv_path = energy_path.string();
    }

    DensityLogConfig density_cfg;
    fs::path density_dir;
    if (P.export_ref_density) {
        density_cfg.enabled = true;
        density_cfg.x_inner = x_inner;
        density_dir = io::make_simulation_subdir(out_dir, P.U0, P.gamma, dt);
        fs::create_directories(density_dir);
        density_cfg.num_csv_path = (density_dir / "num_density.csv").string();
        density_cfg.ref_csv_path = (density_dir / "ref_density.csv").string();
        if (!P.quiet) {
            std::cout << "# [density] exporting numerical/reference densities to "
                      << density_dir << "\n";
        }
    }

    evolve(method, T, spectral, psi_init, g.dx, dt, tol, nsteps, K,
           csv_path.string(), x_ptr,
           P.wide_re, P.wide_im, P.quiet,
           P.log_every, P.csv_every, P.aggregate,
           P.flush_every, P.no_theta, P.profile,
           energy_cfg, density_cfg);

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

