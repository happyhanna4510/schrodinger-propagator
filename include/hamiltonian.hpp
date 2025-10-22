#pragma once
#include <Eigen/Dense>
#include <vector>
#include "grid.hpp"

// Собираем (N-2)x(N-2) матрицу H на внутренних узлах (Дирихле: psi=0 на краях)
Eigen::MatrixXd build_hamiltonian(const Grid& g, const std::vector<double>& U);

// Нормировка собственных векторов по интегралу \sum |phi|^2 dx = 1
void renormalize_eigenvectors(Eigen::MatrixXd& EV, double dx);
