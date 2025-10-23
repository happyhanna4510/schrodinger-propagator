#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

#include "core/io_utils.hpp"
#include "core/spectral.hpp"
#include "core/tridiag.hpp"
#include "core/workspace.hpp"

void taylor_step_tridiag(const Tridiag& T,
                         Eigen::VectorXcd& psi,
                         double dt,
                         int K,
                         TaylorWorkspace& workspace);

void evolve_taylor_tridiag(const Tridiag& T,
                           const SpectralData& spectral,
                           const Eigen::VectorXcd& psi_init,
                           double dx,
                           double dt,
                           int nsteps,
                           int K,
                           int log_every,
                           const std::string& csv_path,
                           const std::vector<double>* x_inner,
                           bool wide_re,
                           bool wide_im,
                           LogExtras extras,
                           bool quiet);

