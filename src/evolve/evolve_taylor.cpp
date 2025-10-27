#include "evolve/evolve_taylor.hpp"

#include <complex>

#include "core/math_utils.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);
}

void taylor_step_tridiag(const Tridiag& T,
                         Eigen::VectorXcd& psi,
                         double dt,
                         int K,
                         TaylorWorkspace& workspace) {
    workspace.resize(psi.size());
    auto& sum = workspace.sum();
    auto& vk  = workspace.vk();
    auto& tmp = workspace.tmp();

    sum = psi;
    vk  = psi;

    const std::complex<double> scale = -I * dt;
    for (int k = 1; k <= K; ++k) {
        tridiag_mul(T, vk, tmp);
        tmp *= (scale / static_cast<double>(k));
        vk = tmp;
        sum.noalias() += vk;
    }

    psi.swap(sum);
}

TaylorEvolver::TaylorEvolver(const Tridiag& T, const EvolverConfig& cfg)
    : EvolverBase(T, cfg) {
    workspace_.resize(T.a.size());
}

StepResult TaylorEvolver::step(Eigen::VectorXcd& psi) {
    workspace_.resize(psi.size());
    taylor_step_tridiag(T_, psi, cfg_.dt, cfg_.K, workspace_);

    StepResult result;
    result.matvecs = cfg_.K;
    return result;
}

