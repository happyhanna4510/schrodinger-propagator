#pragma once
#include <Eigen/Dense>
#include <vector>
#include "grid.hpp"

Eigen::MatrixXd build_hamiltonian(const Grid& g, const std::vector<double>& U);

void renormalize_eigenvectors(Eigen::MatrixXd& EV, double dx);
