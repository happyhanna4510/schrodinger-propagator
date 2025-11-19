#pragma once

#include <Eigen/Dense>

#include <complex>
#include <vector>

#include "core/tridiag.hpp"

std::complex<double> inner_dx(const Eigen::VectorXcd& a,
                              const Eigen::VectorXcd& b,
                              double dx);

double l2_norm_sq(const Eigen::VectorXcd& v, double dx);
double l2_norm(const Eigen::VectorXcd& v, double dx);

double prob_slice(const Eigen::VectorXcd& psi, int i0, int i1, double dx);

double compute_energy(const Eigen::VectorXcd& psi,
                      const std::vector<double>& potential,
                      double dx,
                      double hbar,
                      double mass);

void tridiag_mul(const Tridiag& T,
                 const Eigen::VectorXcd& x,
                 Eigen::VectorXcd& y,
                 int* reallocations = nullptr);

Tridiag make_tridiag_from_dense(const Eigen::Ref<const Eigen::MatrixXd>& H);

void maybe_warn_timestep(double dt, double dE, bool quiet);

