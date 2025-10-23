#pragma once

#include <Eigen/Dense>

struct SpectralData {
    Eigen::VectorXd  evals;
    Eigen::MatrixXcd eigenvectors;
    Eigen::VectorXcd coeffs0;
    Eigen::VectorXcd phi0;
    double           E0 = 0.0;
};

SpectralData make_spectral_data(const Eigen::MatrixXd& H,
                                const Eigen::VectorXcd& psi_init,
                                double dx);

