#include "core/math_utils.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

// tworzy strukturę trójdiagonalną (wektory a,b,c) z gęstej macierzy H
// kopiując diagonalę oraz pod/naddiagonalę
Tridiag make_tridiag_from_dense(const Eigen::Ref<const Eigen::MatrixXd>& H) 
{
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


// oblicza kwadrat normy L2 wektora zespolonego z wagą siatki dx
// sum_i |v_i|^2 * dx 
double l2_norm_sq(const Eigen::VectorXcd& v, double dx) 
{
    long double sum = 0.0L;
    for (Eigen::Index i = 0; i < v.size(); ++i) {
        sum += std::norm(v[i]);
    }
    return static_cast<double>(sum * static_cast<long double>(dx));
}

//pierwiastek z l2_norm_sq
double l2_norm(const Eigen::VectorXcd& v, double dx) {
    return std::sqrt(std::max(0.0, l2_norm_sq(v, dx)));
}

//iloczyn skalarny wektorów zespolonych z wagą dx
std::complex<double> inner_dx(const Eigen::VectorXcd& a,
                              const Eigen::VectorXcd& b,
                              double dx) {
    return dx * a.dot(b);
}

double prob_slice(const Eigen::VectorXcd& psi, int i0, int i1, double dx)
{
    double s = 0.0;
    for (int i = i0; i < i1; ++i) {
        s += std::norm(psi[i]);
    }
    return s * dx;
}

double compute_energy(const Eigen::VectorXcd& psi,
                      const std::vector<double>& potential,
                      double dx,
                      double hbar,
                      double mass)
{
    const int N = static_cast<int>(psi.size());
    if (potential.size() != static_cast<size_t>(N)) {
        throw std::invalid_argument("compute_energy: potential size mismatch with psi");
    }

    std::complex<double> E = 0.0;

    for (int i = 0; i < N; ++i) {
        const std::complex<double> a = (i == 0) ? std::complex<double>(0.0, 0.0) : psi[i - 1];
        const std::complex<double> c = (i == N - 1) ? std::complex<double>(0.0, 0.0) : psi[i + 1];
        const std::complex<double> b = psi[i];

        E += dx * std::conj(b) * (
            -(hbar * hbar) / (2.0 * mass) * (a - 2.0 * b + c) / (dx * dx)
            + potential[static_cast<std::size_t>(i)] * b);
    }

    return E.real();
}

//mnoży macierz trójdiagonalną T przez wektor x i zapisuje wynik w y
void tridiag_mul(const Tridiag& T,
                 const Eigen::VectorXcd& x,
                 Eigen::VectorXcd& y,
                 int* reallocations)
{
    const int M = static_cast<int>(T.a.size());
    if (y.size() != M) {
        if (reallocations) {
            *reallocations += 1;
        }
        y.resize(M);
    }

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

//ostrzeżenie o stabilności kroku czasowego
void maybe_warn_timestep(double dt, double dE, bool quiet) 
{
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

