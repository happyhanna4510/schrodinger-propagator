#include "taylor.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "hamiltonian.hpp"
#include "initial.hpp"
#include "io.hpp"

namespace fs = std::filesystem;

namespace {

constexpr cplx I(0.0, 1.0);

double l2_norm_sq(const Eigen::VectorXcd& v, double dx) {
    long double sum = 0.0L;
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        sum += std::norm(v[i]);
    }
    return static_cast<double>(sum * static_cast<long double>(dx));
}

double l2_norm(const Eigen::VectorXcd& v, double dx) {
    return std::sqrt(std::max(0.0, l2_norm_sq(v, dx)));
}

cplx inner_dx(const Eigen::VectorXcd& a, const Eigen::VectorXcd& b, double dx) {
    return dx * a.dot(b);
}

struct TaylorWorkspace {
    void resize(Eigen::Index n) {
        if (sum_.size() != n) {
            sum_.resize(n);
            vk_.resize(n);
            tmp_.resize(n);
        }
    }

    Eigen::VectorXcd& sum() { return sum_; }
    Eigen::VectorXcd& vk() { return vk_; }
    Eigen::VectorXcd& tmp() { return tmp_; }

private:
    Eigen::VectorXcd sum_;
    Eigen::VectorXcd vk_;
    Eigen::VectorXcd tmp_;
};

void tridiag_mul(const Tridiag& T,
                 const Eigen::VectorXcd& x,
                 Eigen::VectorXcd& y) {
    const int M = static_cast<int>(T.a.size());
    y.resize(M);

    if (M == 0) {
        return;
    }

    if (M == 1) {
        y[0] = T.a[0] * x[0];
        return;
    }

    y[0] = T.a[0] * x[0] + T.b[0] * x[1];
    for (int i = 1; i < M - 1; ++i) {
        y[i] = T.c[i - 1] * x[i - 1] + T.a[i] * x[i] + T.b[i] * x[i + 1];
    }
    y[M - 1] = T.c[M - 2] * x[M - 2] + T.a[M - 1] * x[M - 1];
}

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

    const cplx scale = -I * dt;
    for (int k = 1; k <= K; ++k) {
        tridiag_mul(T, vk, tmp);
        tmp *= (scale / static_cast<double>(k));
        vk = tmp;
        sum.noalias() += vk;
    }

    psi.swap(sum);
}

double prob_slice(const Eigen::VectorXcd& psi, int i0, int i1, double dx) {
    double s = 0.0;
    for (int i = i0; i < i1; ++i) {
        s += std::norm(psi[i]);
    }
    return s * dx;
}

struct SpectralData {
    Eigen::VectorXd  evals;
    Eigen::MatrixXcd eigenvectors; // столбцы — нормированные собственные функции
    Eigen::VectorXcd coeffs0;      // проекция начального состояния
    Eigen::VectorXcd phi0;         // основное состояние
    double           E0 = 0.0;
};

SpectralData make_spectral_data(const Eigen::MatrixXd& H,
                                const Eigen::VectorXcd& psi_init,
                                double dx) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("eigensolve failed");
    }

    Eigen::MatrixXd Vreal = es.eigenvectors();
    renormalize_eigenvectors(Vreal, dx);

    SpectralData data;
    data.evals = es.eigenvalues();
    data.eigenvectors = Vreal.cast<cplx>();

    data.phi0 = data.eigenvectors.col(0);
    if (data.phi0(0).real() < 0.0) {
        data.phi0 = -data.phi0;
        data.eigenvectors.col(0) = -data.eigenvectors.col(0);
    }
    data.E0 = data.evals(0);

    data.coeffs0 = dx * data.eigenvectors.adjoint() * psi_init;

    return data;
}

struct LogSnapshot {
    int    step = 0;
    double t = 0.0;
    double norm_sq = 0.0;
    double norm = 0.0;
    double drift_sq = 0.0;
    double drift = 0.0;
    double p0 = 0.0;
    double err_exact = 0.0;
    double err_true = 0.0;
    double theta = 0.0;
    double edge_left = 0.0;
    double edge_right = 0.0;
};

LogSnapshot collect_snapshot(const SpectralData& spectral,
                             const Eigen::VectorXcd& psi,
                             double t,
                             double dx,
                             int step,
                             double norm_sq0,
                             double norm0) {
    LogSnapshot snap;
    snap.step = step;
    snap.t = t;
    snap.norm_sq = l2_norm_sq(psi, dx);
    snap.norm = std::sqrt(std::max(0.0, snap.norm_sq));
    snap.drift_sq = snap.norm_sq - norm_sq0;
    snap.drift = snap.norm - norm0;

    snap.p0 = std::norm(inner_dx(spectral.phi0, psi, dx));

    const cplx phase0 = std::exp(-I * spectral.E0 * t);
    Eigen::VectorXcd phi_exact = phase0 * spectral.phi0;
    snap.err_exact = (psi - phi_exact).norm();

    Eigen::ArrayXcd phase = (-I * spectral.evals.array() * t).exp();
    Eigen::VectorXcd coeff_exact = (spectral.coeffs0.array() * phase).matrix();
    Eigen::VectorXcd psi_exact = spectral.eigenvectors * coeff_exact;
    snap.err_true = (psi - psi_exact).norm();

    Eigen::VectorXcd coeff_num = dx * spectral.eigenvectors.adjoint() * psi;
    snap.theta = (coeff_exact - coeff_num).squaredNorm();

    const int width = std::max(1, std::min<int>(10, static_cast<int>(psi.size())));
    snap.edge_left = prob_slice(psi, 0, width, dx);
    snap.edge_right = prob_slice(psi, static_cast<int>(psi.size()) - width,
                                 static_cast<int>(psi.size()), dx);

    return snap;
}

void write_csv_header(std::ofstream& f, LogExtras extras) {
    f << "step,t,norm,norm_sqrt";
    if (has(extras, LogExtras::P0)) {
        f << ",p0";
    }
    if (has(extras, LogExtras::ErrPhi)) {
        f << ",err_exact";
    }
    f << ",err_true,theta,Pleft,Pright\n";
}

void write_csv_row(std::ofstream& f, const LogSnapshot& snap, LogExtras extras) {
    std::ostringstream row;
    row << snap.step << ','
        << std::setprecision(15) << std::defaultfloat << snap.t << ','
        << std::setprecision(15) << snap.norm_sq << ','
        << std::setprecision(15) << snap.norm;

    if (has(extras, LogExtras::P0)) {
        row << ',' << std::setprecision(12) << snap.p0;
    }
    if (has(extras, LogExtras::ErrPhi)) {
        row << ',' << std::scientific << std::setprecision(6) << snap.err_exact
            << std::defaultfloat;
    }

    row << ',' << std::scientific << std::setprecision(6) << snap.err_true
        << ',' << std::scientific << std::setprecision(6) << snap.theta
        << ',' << std::scientific << std::setprecision(6) << snap.edge_left
        << ',' << std::scientific << std::setprecision(6) << snap.edge_right
        << '\n';

    f << row.str();
}

void print_snapshot(const LogSnapshot& snap,
                    LogExtras extras) {
    std::ostringstream line;
    line << std::setw(8) << snap.step
         << "  t=" << std::scientific << std::setprecision(6) << snap.t
         << std::defaultfloat
         << "  norm=" << std::setprecision(15) << snap.norm_sq
         << "  |psi|=" << std::fixed << std::setprecision(12) << snap.norm
         << std::scientific
         << "  drift=" << std::setprecision(3) << snap.drift_sq
         << "  drift|psi|=" << std::setprecision(3) << snap.drift;

    line << std::defaultfloat;

    if (has(extras, LogExtras::P0)) {
        line << "  p0=" << std::fixed << std::setprecision(10) << snap.p0;
    }
    if (has(extras, LogExtras::ErrPhi)) {
        line << std::scientific << "  err_exact=" << std::setprecision(3) << snap.err_exact;
    }

    line << std::scientific
         << "  err_true=" << std::setprecision(3) << snap.err_true
         << "  Theta=" << std::setprecision(3) << snap.theta
         << '\n'
         << "  edge(P): L=" << std::setprecision(3) << snap.edge_left
         << " R=" << std::setprecision(3) << snap.edge_right;

    std::cout << line.str() << '\n';
}

void maybe_warn_timestep(double dt, double dE, bool quiet) {
    if (quiet) {
        return;
    }
    if (dt <= 0.0 || dE <= 0.0) {
        return;
    }

    const double stability = dE * dt;
    if (stability > 0.5) {
        std::cout << "# [warn] dt*dE ≈ " << stability
                  << ": consider reducing dt or increasing K for stability\n";
    }
}

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

} // namespace

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

fs::path run_taylor_evolution(const Grid& g,
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

    fs::create_directories(out_dir);
    fs::path csv_path = P.csv_name.empty()
        ? io::make_csv_path(out_dir, P, K, dt, g)
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

    if (!P.quiet) {
        std::cout << "\n# === Time evolution (Taylor K=" << K
                  << ", dt=" << dt << ", tmax=" << tmax << ") ===\n";
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

    evolve_taylor_tridiag(T, spectral, psi_init, g.dx, dt, nsteps, K,
                          log_every, csv_path.string(), x_ptr,
                          P.wide_re, P.wide_im, extras, P.quiet);

    if (!P.quiet) {
        std::cout << "# log saved to: " << csv_path.string() << "\n";
    }

    return csv_path;
}

