#include "evolve/evolve_rk4.hpp"

#include <complex>

#include "core/math_utils.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);
}

Rk4Evolver::Rk4Evolver(const Tridiag& T, const EvolverConfig& cfg)
    : EvolverBase(T, cfg) 
    {
    const Eigen::Index n = T.a.size();
    k1_.resize(n);
    k2_.resize(n);
    k3_.resize(n);
    k4_.resize(n);
    tmp_.resize(n);
}

StepResult Rk4Evolver::step(Eigen::VectorXcd& psi) 
{
    if (k1_.size() != psi.size()) {
        const Eigen::Index n = psi.size();
        k1_.resize(n);
        k2_.resize(n);
        k3_.resize(n);
        k4_.resize(n);
        tmp_.resize(n);
    }

    tridiag_mul(T_, psi, k1_);
    k1_ *= -I;

    tmp_ = psi + (0.5 * cfg_.dt) * k1_;
    tridiag_mul(T_, tmp_, k2_);
    k2_ *= -I;

    tmp_ = psi + (0.5 * cfg_.dt) * k2_;
    tridiag_mul(T_, tmp_, k3_);
    k3_ *= -I;

    tmp_ = psi + cfg_.dt * k3_;
    tridiag_mul(T_, tmp_, k4_);
    k4_ *= -I;

    Eigen::VectorXcd sum = k1_;
    sum.noalias() += 2.0 * k2_;
    sum.noalias() += 2.0 * k3_;
    sum.noalias() += k4_;
    psi.noalias() += (cfg_.dt / 6.0) * sum;

    StepResult result;
    result.matvecs = 4; //rk4 zawsze używa 4 mnożeń macierz–wektor
    return result;
}

