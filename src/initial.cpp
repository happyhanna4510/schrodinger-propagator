#include "initial.hpp"
#include <cmath>
#include <complex>

// znormalizowana (L2 z wagą dx) gaussowska paczka na punktach wewnętrznych siatki
// x0 środek, sigma szerokość
Eigen::VectorXd gaussian_on_inner(const Grid& g, double x0, double sigma)
{
    const Eigen::Index M = static_cast<Eigen::Index>(g.N - 2);  // liczba punktów wewnętrznych
    Eigen::VectorXd v(M);
    for (Eigen::Index i = 0; i < M; ++i) {
        double x = g.x[static_cast<Eigen::Index>(i + 1)]; // pomijamy punkty brzegowe
        v(i) = std::exp(-(x - x0) * (x - x0) / (2 * sigma * sigma));
    }
    double norm = (v.array().square().sum()) * g.dx;
    v /= std::sqrt(norm);
    return v;
}

// (zespolona) gaussowska paczka z fazą fali płaskiej e^{i k0 x}
// wygodna do ruchomej paczki (pęd ~ k0).
using cplx = std::complex<double>;
Eigen::VectorXcd gaussian_complex_on_inner(const Grid& g, double x0, double sigma, double k0)
{
    Eigen::VectorXcd v = gaussian_on_inner(g, x0, sigma).cast<cplx>();
    const Eigen::Index M = v.size();

    for (Eigen::Index i = 0; i < M; ++i) {
        double x = g.x[static_cast<Eigen::Index>(i + 1)];
        v(i) *= std::exp(cplx(0.0, k0 * x));
    }

    double norm = 0.0;
    for (Eigen::Index i = 0; i < M; ++i) {
        norm += std::norm(v(i));
    }
    v /= std::sqrt(norm * g.dx);
    return v;
}
