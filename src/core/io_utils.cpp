#include "core/io_utils.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <cstdlib>
#include <utility>

#include "core/math_utils.hpp"

namespace 
{
    //i = sqrt(-1)
constexpr std::complex<double> I(0.0, 1.0);

//formatuje liczbę mnożeń macierz–wektor (matvecs)
std::string format_matvecs(double value) 
{
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

// SFINAE: próbuje pobrać S.dx jeśli istnieje
template <typename Spectral>
auto try_get_dx(const Spectral& S, int) -> decltype(S.dx, double()) {
    return S.dx;
}

// Fallback: gdy S.dx nie istnieje, zwraca NaN.
template <typename Spectral>
double try_get_dx(const Spectral&, ...) {
    return std::numeric_limits<double>::quiet_NaN();
}
}

// relative=false: zwracaj błąd absolutny (||C_dok - C_num||^2)
//relative=true: zwracaj błąd względny (sqrt(theta_abs / ||C_dok||^2)
double compute_theta(const SpectralData& S,
                     const Eigen::VectorXcd& psi,
                     double t, double dx,
                     bool relative /*= false*/) 
{
    static const std::complex<double> iC(0.0, 1.0);

    const bool size_match =
        (S.eigenvectors.rows() == psi.size()) &&
        (S.eigenvectors.cols() == S.evals.size()) &&
        (S.eigenvectors.cols() == S.coeffs0.size());
    if (!size_match) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //S.dx (jeśli jest i >0)->  dx -> 1.0
    auto select_dx = [&](double provided_dx) 
    {
        double dx_eff = provided_dx;

        const double candidate = try_get_dx(S, std::numeric_limits<double>::quiet_NaN());

        if (std::isfinite(candidate) && candidate > 0.0) {
            dx_eff = candidate;
        }
        if (!std::isfinite(dx_eff) || dx_eff <= 0.0) {
            dx_eff = 1.0;
        }
        return dx_eff;
    };
    const double dx_eff = select_dx(dx);

    const Eigen::ArrayXcd phase = (-iC * S.evals.array() * t).exp();
    const Eigen::VectorXcd C_dok = (S.coeffs0.array() * phase).matrix();
    const Eigen::VectorXcd C_num = dx_eff * S.eigenvectors.adjoint() * psi;

    const double theta_abs = (C_dok - C_num).squaredNorm();

    if (!relative) {
        return theta_abs;
    } else {
        const double denom = std::max(1e-300, C_dok.squaredNorm());
        return std::sqrt(theta_abs / denom); 
    }
}



void write_step_csv_header(std::ofstream& f, bool include_cheb_extras) 
{
    f << "method,step,t,dt,dt_ms,matvecs,norm_err,theta_abs,theta_rel";
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
                        double norm_err,
                        double theta_rel,  
                        double theta_abs,  
                        std::optional<int> K_used,
                        std::optional<double> bn_ratio,
                        bool include_cheb_extras) 
{
    std::ostringstream row;
    row.setf(std::ios::fmtflags(0), std::ios::floatfield);

    row << method << ','
        << step << ','
        << std::fixed << std::setprecision(6) << t << ','
        << std::setprecision(15) << std::defaultfloat << dt << ','
        << std::fixed << std::setprecision(3) << dt_ms << ','
        << std::defaultfloat << matvecs << ','

        << std::scientific << std::setprecision(6) << norm_err << ','
        << std::scientific << std::setprecision(6) << theta_rel << ',' 
        << std::scientific << std::setprecision(6) << theta_abs;     

    if (include_cheb_extras) {
        row << ',';
        if (K_used) row << *K_used;
        row << ',';
        if (bn_ratio) row << std::scientific << std::setprecision(6) << *bn_ratio;
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
                        double theta_rel,
                        double theta_abs,
                        std::optional<int> K_used) 
{
    std::ostringstream line;
    line << '[' << method << "] "
         << "step=" << step
         << " t=" << std::fixed << std::setprecision(6) << t
         << "  dt_ms=" << std::fixed << std::setprecision(3) << dt_ms
         << "  matvecs=" << format_matvecs(matvecs)
         << "  norm_err=" << std::scientific << std::setprecision(3) << norm_err
         << "  theta_rel=" << std::scientific << std::setprecision(3) << theta_abs
         << "  theta_abs=" << std::scientific << std::setprecision(3) << theta_rel;


    if (K_used) {
        line << "  K=" << *K_used;
    }

    std::cout << line.str() << '\n';
}

