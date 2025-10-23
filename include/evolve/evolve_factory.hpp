#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

#include "core/io_utils.hpp"
#include "core/spectral.hpp"
#include "core/tridiag.hpp"

void evolve(const std::string& method,
            const Tridiag& T,
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

