#include "core/io_utils.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "core/math_utils.hpp"
#include "core/spectral.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);
}

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

    const std::complex<double> phase0 = std::exp(-I * spectral.E0 * t);
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

