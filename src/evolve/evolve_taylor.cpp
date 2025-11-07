#include "evolve/evolve_taylor.hpp"

#include <chrono>
#include <complex>

#include "core/math_utils.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);
}

void taylor_step_tridiag(const Tridiag& T,
                         Eigen::VectorXcd& psi,
                         double dt,
                         int K,
                         TaylorWorkspace& workspace,
                         StepProfile* profile)
{
    const bool resized = workspace.resize(psi.size());
    auto& sum = workspace.sum();
    auto& vk  = workspace.vk();
    auto& tmp = workspace.tmp();

    sum = psi;
    vk  = psi;

    const std::complex<double> scale = -I * dt;
    std::chrono::steady_clock::time_point series_start;
    if (profile) {
        series_start = std::chrono::steady_clock::now();
        profile->reallocations += resized ? 3 : 0;
    }
    std::chrono::duration<double, std::micro> rhs_accum(0.0);

    for (int k = 1; k <= K; ++k) {
        std::chrono::steady_clock::time_point rhs_start;
        if (profile) {
            rhs_start = std::chrono::steady_clock::now();
        }
        tridiag_mul(T, vk, tmp, profile ? &profile->reallocations : nullptr);
        if (profile) {
            rhs_accum += std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - rhs_start);
        }
        tmp *= (scale / static_cast<double>(k));
        sum.noalias() += tmp;
        vk.swap(tmp);
    }

    psi.swap(sum);

    if (profile) {
        const auto series_end = std::chrono::steady_clock::now();
        profile->work_us = std::chrono::duration<double, std::micro>(series_end - series_start).count();
        profile->rhs_us = rhs_accum.count();
    }
}

TaylorEvolver::TaylorEvolver(const Tridiag& T, const EvolverConfig& cfg)
    : EvolverBase(T, cfg) {
    workspace_.resize(T.a.size());
}

StepResult TaylorEvolver::step(Eigen::VectorXcd& psi)
 {
    StepResult result;
    auto profile = cfg_.profile ? std::make_optional<StepProfile>() : std::nullopt;
    std::chrono::steady_clock::time_point step_start;
    if (profile) {
        step_start = std::chrono::steady_clock::now();
    }

    taylor_step_tridiag(T_, psi, cfg_.dt, cfg_.K, workspace_, profile ? &(*profile) : nullptr);

    if (profile) {
        profile->step_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - step_start).count();
        result.profile = std::move(profile);
    }

    result.matvecs = cfg_.K;
    return result;
}

