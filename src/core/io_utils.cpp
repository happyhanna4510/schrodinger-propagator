#include "core/io_utils.hpp"

#include <algorithm>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "core/math_utils.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);

std::string format_matvecs(double value) {
    std::ostringstream oss;
    if (!std::isfinite(value)) {
        oss << "nan";
        return oss.str();
    }
    const double rounded = std::round(value);
    if (std::abs(value - rounded) < 1e-9) {
        oss << static_cast<long long>(std::llround(rounded));
    } else {
        oss << std::fixed << std::setprecision(3) << value;
    }
    return oss.str();
}
}

double compute_theta(const SpectralData& S,
                     const Eigen::VectorXcd& psi,
                     double t, double dx,
                     bool relative /*= false*/) {
    using cplx = std::complex<double>;

    if (!std::isfinite(dx)) std::cerr << "[θ] dx=" << dx << " (bad)\n";
    if (!S.eigenvectors.allFinite()) std::cerr << "[θ] eigenvectors NaN\n";
    if (!S.coeffs0.allFinite()) std::cerr << "[θ] coeffs0 NaN\n";
    if (!S.evals.allFinite()) std::cerr << "[θ] evals NaN\n";

    std::cerr << "[θ] psi=" << psi.size()
          << " eigV rows=" << S.eigenvectors.rows()
          << " cols=" << S.eigenvectors.cols()
          << " coeffs0=" << S.coeffs0.size()
          << " evals=" << S.evals.size() << "\n";


    static const cplx iC(0.0, 1.0);

    // --- страховка dx ---
    if (!std::isfinite(dx) || dx <= 0.0) {
        // если у тебя есть S.dx — лучше взять его; иначе хотя бы 1.0
        // dx = S.dx; // если поле есть
        dx = 1.0;
        // временный вывод, чтобы увидеть проблему у источника:
        std::cerr << "[theta] bad dx, fallback to " << dx << "\n";
    }

    // точные коэффициенты: C_dok(t) = C(0) * exp(-i E t)
    Eigen::ArrayXcd phase = (-iC * S.evals.array() * t).exp();
    Eigen::VectorXcd C_dok = (S.coeffs0.array() * phase).matrix();

    // численные: C_num(t) = dx * V^H * ψ(t)
    Eigen::VectorXcd C_num = dx * S.eigenvectors.adjoint() * psi;

    double num = (C_dok - C_num).squaredNorm();
    if (!relative) return num;

    double den = std::max(1e-300, C_dok.squaredNorm());

    if (relative) {
    std::cerr << "[θ] num=" << num
              << " den=" << C_dok.squaredNorm()
              << " rel=" << std::sqrt(num / std::max(1e-300, C_dok.squaredNorm()))
              << "\n";
}

    return std::sqrt(num / den);
}


std::optional<double> compute_e_true(const SpectralData& spectral,
                                     const Eigen::VectorXcd& psi,
                                     double t) {
    if (spectral.evals.size() == 0 || spectral.coeffs0.size() == 0 ||
        spectral.eigenvectors.size() == 0) {
        return std::nullopt;
    }

    const auto n_states = spectral.eigenvectors.cols();
    if (spectral.eigenvectors.rows() != psi.size() ||
        spectral.coeffs0.size() != n_states ||
        spectral.evals.size() != n_states) {
        return std::nullopt;
    }

    Eigen::ArrayXcd phase = (-I * spectral.evals.array() * t).exp();
    Eigen::VectorXcd coeff_exact = (spectral.coeffs0.array() * phase).matrix();
    Eigen::VectorXcd psi_exact = spectral.eigenvectors * coeff_exact;

    const double denom = psi_exact.norm();
    const double numer = (psi - psi_exact).norm();
    if (denom > 0.0) {
        return numer / denom;
    }
    return numer;
}

void write_step_csv_header(std::ofstream& f, bool include_cheb_extras) {
    f << "method,step,t,dt,dt_ms,matvecs,norm_err,theta";
    if (include_cheb_extras) {
        f << ",K_used,bn_ratio";
    }
    f << ",e_true\n";
}

void write_step_csv_row(std::ofstream& f,
                        const std::string& method,
                        int step,
                        double t,
                        double dt,
                        double dt_ms,
                        double matvecs,
                        double norm_err,
                        double theta,
                        std::optional<double> e_true,
                        std::optional<int> K_used,
                        std::optional<double> bn_ratio,
                        bool include_cheb_extras) {
    std::ostringstream row;
    row << method << ','
        << step << ','
        << std::fixed << std::setprecision(6) << t << ','
        << std::setprecision(15) << std::defaultfloat << dt << ','
        << std::fixed << std::setprecision(3) << dt_ms << ','
        << std::defaultfloat << matvecs << ','
        << std::scientific << std::setprecision(6) << norm_err << ','
        << std::scientific << std::setprecision(6) << theta;

    if (include_cheb_extras) {
        if (K_used) {
            row << ',' << *K_used;
        } else {
            row << ',';
        }
        if (bn_ratio) {
            row << ',' << std::scientific << std::setprecision(6) << *bn_ratio;
        } else {
            row << ',';
        }
    }

    row << ',';
    if (e_true) {
        row << std::scientific << std::setprecision(6) << *e_true;
    }

    row << '\n';
    f << row.str();
}

void print_step_console(const std::string& method,
                        int step,
                        double t,
                        double dt_ms,
                        double matvecs,
                        double norm_err,
                        double theta,
                        std::optional<double> e_true,
                        std::optional<int> K_used) {
    std::ostringstream line;
    line << '[' << method << "] "
         << "step=" << step
         << " t=" << std::fixed << std::setprecision(6) << t
         << "  dt_ms=" << std::fixed << std::setprecision(3) << dt_ms
         << "  matvecs=" << format_matvecs(matvecs)
         << "  norm_err=" << std::scientific << std::setprecision(3) << norm_err
         << "  theta=" << std::scientific << std::setprecision(3) << theta;

    if (e_true) {
        line << "  Etrue=" << std::scientific << std::setprecision(3) << *e_true;
    }

    if (K_used) {
        line << "  K=" << *K_used;
    }

    std::cout << line.str() << '\n';
}

