#pragma once

#include <Eigen/Core>

#include "core/tridiag.hpp"

struct RevTestResult {
    double err_back = 0.0;
    double norm0 = 0.0;
    double norm_fwd = 0.0;
    double norm_back = 0.0;
    double max_drift = 0.0;
};

RevTestResult test_time_reversibility(const Tridiag& T,
                                      const Eigen::VectorXcd& psi_init,
                                      double dx,
                                      double dt,
                                      int nsteps,
                                      int K);

