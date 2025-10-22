#pragma once
#include "grid.hpp"
#include <Eigen/Dense>

Eigen::VectorXd gaussian_on_inner(const Grid& g, double x0, double sigma);

Eigen::VectorXcd gaussian_complex_on_inner(const Grid& g, double x0, double sigma, double k0);