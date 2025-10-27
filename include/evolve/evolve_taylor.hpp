#pragma once

#include <Eigen/Core>

#include <string>

#include "core/workspace.hpp"
#include "evolve/evolver_base.hpp"

void taylor_step_tridiag(const Tridiag& T,
                         Eigen::VectorXcd& psi,
                         double dt,
                         int K,
                         TaylorWorkspace& workspace);

class TaylorEvolver : public EvolverBase {
public:
    TaylorEvolver(const Tridiag& T, const EvolverConfig& cfg);

    StepResult step(Eigen::VectorXcd& psi) override;

private:
    TaylorWorkspace workspace_;
};

