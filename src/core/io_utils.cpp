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
    const double rounded = std::round(value);
    if (std::abs(value - rounded) < 1e-9) {
        oss << static_cast<long long>(std::llround(rounded));
    } else {
        oss << std::fixed << std::setprecision(3) << value;
    }
    return oss.str();
}
}

StepMetrics compute_step_metrics(const SpectralData& spectral,
                                 const Eigen::VectorXcd& psi,
                                 double t,
                                 double dx) {
    StepMetrics metrics;
    metrics.norm = l2_norm_sq(psi, dx);
    metrics.norm_err = std::abs(metrics.norm - 1.0);

    Eigen::ArrayXcd phase = (-I * spectral.evals.array() * t).exp();
    Eigen::VectorXcd coeff_exact = (spectral.coeffs0.array() * phase).matrix();
    Eigen::VectorXcd coeff_num = dx * spectral.eigenvectors.adjoint() * psi;
    metrics.theta = (coeff_exact - coeff_num).squaredNorm();

    return metrics;
}

void write_step_csv_header(std::ofstream& f, bool include_cheb_extras) {
    f << "method,step,t,dt,dt_ms,matvecs,norm,norm_err,theta";
    if (include_cheb_extras) {
        f << ",K_used,bn_ratio";
    }
    f << '\n';
}

void write_step_csv_row(std::ofstream& f,
                        const std::string& method,
                        int step,
                        double t,
                        double dt,
                        double dt_ms,
                        double matvecs,
                        const StepMetrics& metrics,
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
        << std::setprecision(15) << metrics.norm << ','
        << std::scientific << std::setprecision(6) << metrics.norm_err << ','
        << std::scientific << std::setprecision(6) << metrics.theta;

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

    row << '\n';
    f << row.str();
}

void print_step_console(const std::string& method,
                        int step,
                        double t,
                        double dt_ms,
                        double matvecs,
                        const StepMetrics& metrics,
                        std::optional<int> K_used) {
    std::ostringstream line;
    line << '[' << method << "] "
         << "step=" << step
         << " t=" << std::fixed << std::setprecision(6) << t
         << "  dt_ms=" << std::fixed << std::setprecision(3) << dt_ms
         << "  matvecs=" << format_matvecs(matvecs)
         << "  norm_err=" << std::scientific << std::setprecision(3) << metrics.norm_err
         << "  theta=" << std::scientific << std::setprecision(3) << metrics.theta;

    if (K_used) {
        line << "  K=" << *K_used;
    }

    std::cout << line.str() << '\n';
}

