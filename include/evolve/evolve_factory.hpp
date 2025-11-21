#pragma once

#include <Eigen/Core>

#include <string>
#include <vector>

#include "core/spectral.hpp"
#include "core/tridiag.hpp"

struct EnergyLogConfig {
    bool enabled = false;
    const std::vector<double>* potential = nullptr;
    double hbar = 1.0;
    double mass = 1.0;
    std::string csv_path;
};

struct DensityLogConfig {
    bool enabled = false;
    std::vector<double> x_inner;
    std::string num_csv_path;
    std::string ref_csv_path;
};

void evolve(const std::string& method,
            const Tridiag& T,
            const SpectralData& spectral,
            const Eigen::VectorXcd& psi_init,
            double dx,
            double dt,
            double tol,
            int nsteps,
            int K,
            const std::string& csv_path,
            const std::vector<double>* x_inner,
            bool wide_re,
            bool wide_im,
            bool quiet,
            int log_every,
            int csv_every,
            bool aggregate,
            int flush_every,
            bool no_theta,
            bool profile,
            const EnergyLogConfig& energy_cfg = EnergyLogConfig{},
            const DensityLogConfig& density_cfg = DensityLogConfig{});

