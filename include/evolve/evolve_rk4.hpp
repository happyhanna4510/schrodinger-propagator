#pragma once

#include <Eigen/Core>

#include "evolve/evolver_base.hpp"

class Rk4Evolver : public EvolverBase {
public:
    Rk4Evolver(const Tridiag& T, const EvolverConfig& cfg);

    StepResult step(Eigen::VectorXcd& psi) override;

private:
    Eigen::VectorXcd k1_;
    Eigen::VectorXcd k2_;
    Eigen::VectorXcd k3_;
    Eigen::VectorXcd k4_;
    Eigen::VectorXcd tmp_;
};

