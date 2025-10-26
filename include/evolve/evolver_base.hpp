#pragma once

#include <Eigen/Core>

#include <optional>

#include "core/tridiag.hpp"

struct EvolverConfig {
    double dt        = 0.0;
    double dx        = 0.0;
    double hbar      = 1.0;
    double tolerance = 1e-12;
    int    K         = 0;
    int    log_every = 10000;
    int    csv_every = 1;
    bool   aggregate = true;
};

struct StepResult {
    int matvecs = 0;
    std::optional<int>    K_used;
    std::optional<double> bn_ratio;
};

class EvolverBase {
public:
    EvolverBase(const Tridiag& T, const EvolverConfig& cfg)
        : T_(T)
        , cfg_(cfg) {}

    virtual ~EvolverBase() = default;

    virtual StepResult step(Eigen::VectorXcd& psi) = 0;

protected:
    const Tridiag&  T_;
    EvolverConfig   cfg_;
};

