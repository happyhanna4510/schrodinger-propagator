#include "evolve/evolve_rk4.hpp"

#include <chrono>
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
    sum_.resize(n);
}

StepResult Rk4Evolver::step(Eigen::VectorXcd& psi)
{
    StepResult result;
    auto profile = cfg_.profile ? std::make_optional<StepProfile>() : std::nullopt;
    std::chrono::steady_clock::time_point step_start;
    if (profile) {
        step_start = std::chrono::steady_clock::now();
    }

    if (k1_.size() != psi.size()) {
        const Eigen::Index n = psi.size();
        k1_.resize(n);
        k2_.resize(n);
        k3_.resize(n);
        k4_.resize(n);
        tmp_.resize(n);
        sum_.resize(n);
        if (profile) {
            profile->reallocations += 6;
        }
    }

    std::chrono::duration<double, std::micro> rhs_accum(0.0);

    std::chrono::steady_clock::time_point rhs_start;
    if (profile) rhs_start = std::chrono::steady_clock::now();
    tridiag_mul(T_, psi, k1_, profile ? &profile->reallocations : nullptr);
    if (profile) rhs_accum += std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - rhs_start);
    k1_ *= -I;

    const double dt = cfg_.dt;

    tmp_ = psi;
    tmp_.noalias() += (0.5 * dt) * k1_;
    if (profile) rhs_start = std::chrono::steady_clock::now();
    tridiag_mul(T_, tmp_, k2_, profile ? &profile->reallocations : nullptr);
    if (profile) rhs_accum += std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - rhs_start);
    k2_ *= -I;

    tmp_ = psi;
    tmp_.noalias() += (0.5 * dt) * k2_;
    if (profile) rhs_start = std::chrono::steady_clock::now();
    tridiag_mul(T_, tmp_, k3_, profile ? &profile->reallocations : nullptr);
    if (profile) rhs_accum += std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - rhs_start);
    k3_ *= -I;

    tmp_ = psi;
    tmp_.noalias() += dt * k3_;
    if (profile) rhs_start = std::chrono::steady_clock::now();
    tridiag_mul(T_, tmp_, k4_, profile ? &profile->reallocations : nullptr);
    if (profile) rhs_accum += std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - rhs_start);
    k4_ *= -I;

    sum_.noalias() = k1_;
    sum_.noalias() += (2.0 * k2_);
    sum_.noalias() += (2.0 * k3_);
    sum_.noalias() += k4_;
    psi.noalias() += (dt / 6.0) * sum_;

    result.matvecs = 4; //rk4 zawsze używa 4 mnożeń macierz–wektor
    if (profile) {
        profile->rhs_us = rhs_accum.count();
        profile->step_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - step_start).count();
        profile->work_us = profile->step_us - profile->rhs_us;
        result.profile = std::move(profile);
    }
    return result;
}

