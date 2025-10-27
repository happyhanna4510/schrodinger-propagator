#pragma once

#include <Eigen/Core>

#include <complex>

#include "core/spectral.hpp"
#include "evolve/evolver_base.hpp"

class ChebyshevEvolver : public EvolverBase {
public:
    ChebyshevEvolver(const Tridiag& T,
                     const SpectralData& spectral,
                     const EvolverConfig& cfg);

    StepResult step(Eigen::VectorXcd& psi) override;

private:
    double Emin_ = 0.0;
    double Emax_ = 0.0;
    double deltaE_ = 0.0;
    double center_ = 0.0;
    double norm_factor_ = 0.0;
    double x_ = 0.0;
    std::complex<double> phase_ = {1.0, 0.0};
    std::complex<double> trivial_phase_ = {1.0, 0.0};
    bool   trivial_case_ = false;

    Eigen::VectorXcd p_prev_;
    Eigen::VectorXcd p_curr_;
    Eigen::VectorXcd tmp_;

    void apply_normalized(const Eigen::VectorXcd& in, Eigen::VectorXcd& out);
};

