#include "morse_static.hpp"
#include "hamiltonian.hpp"
#include "morse_analytic.hpp"
#include "io.hpp"
#include <Eigen/Eigenvalues>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

void write_morse_outputs(const Grid& g,
                                const std::vector<double>& U_true,
                                double gamma, int first,
                                const fs::path& out,
                                bool quiet)
{
    std::error_code ec;
    fs::create_directories(out, ec);

    Eigen::MatrixXd H_true = build_hamiltonian(g, U_true);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_true(H_true);
    if (es_true.info()!=Eigen::Success) { std::cerr << "eigensolve failed\n"; return; }

    Eigen::VectorXd Evals_true = es_true.eigenvalues();
    Eigen::MatrixXd EV_true    = es_true.eigenvectors();
    renormalize_eigenvectors(EV_true, g.dx);

    if (!quiet) {
        std::cout << "# Grid: N="<<g.N<<" xmax="<<g.xmax<<" dx="<<g.dx<<"\n";
        std::cout << "# Morse gamma="<<gamma<<"  (n_max="<<morse_nmax(gamma)<<")\n\n";
    }

    const Eigen::Index available = Evals_true.size();
    const Eigen::Index show = std::min<Eigen::Index>(static_cast<Eigen::Index>(first), available);
    auto Eanal = morse_energies(gamma, static_cast<int>(show));
    std::vector<double> Eval_num(static_cast<std::size_t>(show));
    for (Eigen::Index i = 0; i < show; ++i) {
        Eval_num[static_cast<std::size_t>(i)] = Evals_true(i);
    }
    if (!quiet) {
        print_energy_table(Eval_num, Eanal, static_cast<int>(show));
    }

    Eigen::MatrixXd G = (EV_true.transpose() * (g.dx * EV_true));
    double max_dev = (G - Eigen::MatrixXd::Identity(G.rows(), G.cols())).cwiseAbs().maxCoeff();
    if (!quiet) {
        std::cout << "# Orthonorm check (max |G-I|): " << max_dev << "\n";
    }

    std::vector<double> xv(g.N), Uv(g.N);
    for (int i=0;i<g.N;i++){ xv[i]=g.x[i]; Uv[i]=U_true[i]; }
    save_xy_csv(out / "morse_potential.csv", xv, Uv, "x", "U");

    int Ksave = std::min(first, g.N-2);
    std::vector<double> xin(g.N-2);
    for (int i=0;i<g.N-2;i++) xin[i]=g.x[i+1];
    save_matrix_csv(out / "morse_eigenstates.csv", xin, EV_true, Ksave, "x", "psi_n");

    std::vector<double> Eall(static_cast<std::size_t>(Evals_true.size()));
    for (Eigen::Index i = 0; i < Evals_true.size(); ++i) {
        Eall[static_cast<std::size_t>(i)] = Evals_true(i);
    }
    save_vector_csv(out / "morse_energies_num.csv",  Eall, "E_num");

    auto Eanal_all = morse_energies(gamma, morse_nmax(gamma)+1);
    save_vector_csv(out / "morse_energies_anal.csv", Eanal_all, "E_anal");

    if (!quiet) {
        std::cout << "# Morse: morse_* files (re)written in " << out.string() << "\n";
    }
}
