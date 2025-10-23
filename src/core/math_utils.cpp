#include "math_utils.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

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

std::complex<double> inner_dx(const Eigen::VectorXcd& a,
                              const Eigen::VectorXcd& b,
                              double dx) {
    return dx * a.dot(b);
}

double prob_slice(const Eigen::VectorXcd& psi, int i0, int i1, double dx) {
    double s = 0.0;
    for (int i = i0; i < i1; ++i) {
        s += std::norm(psi[i]);
    }
    return s * dx;
}

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
