#pragma once

#include <Eigen/Core>

#include <complex>

#include "taylor.hpp"

double l2_norm_sq(const Eigen::VectorXcd& v, double dx);
double l2_norm(const Eigen::VectorXcd& v, double dx);

std::complex<double> inner_dx(const Eigen::VectorXcd& a,
                              const Eigen::VectorXcd& b,
                              double dx);

double prob_slice(const Eigen::VectorXcd& psi, int i0, int i1, double dx);

void tridiag_mul(const Tridiag& T,
                 const Eigen::VectorXcd& x,
                 Eigen::VectorXcd& y);

void maybe_warn_timestep(double dt, double dE, bool quiet);
