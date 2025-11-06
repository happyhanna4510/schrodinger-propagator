#pragma once
#include <Eigen/Dense>
#include <vector>
#include "grid.hpp"

Eigen::MatrixXd build_hamiltonian(const Grid& g, const std::vector<double>& U, double U0 = 0.0);

void renormalize_eigenvectors(Eigen::MatrixXd& EV, double dx);
