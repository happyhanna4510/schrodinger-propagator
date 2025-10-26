#include "evolve/evolve_chebyshev.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

#include "core/math_utils.hpp"

namespace {
constexpr std::complex<double> I(0.0, 1.0);
constexpr double EPS = 1e-14;
}

ChebyshevEvolver::ChebyshevEvolver(const Tridiag& T,
                                   const SpectralData& spectral,
                                   const EvolverConfig& cfg)
    : EvolverBase(T, cfg) {
    Emin_ = spectral.evals.minCoeff();
    Emax_ = spectral.evals.maxCoeff();
    deltaE_ = Emax_ - Emin_;
    center_ = 0.5 * (Emax_ + Emin_);

    const double scale_ref = std::max({1.0, std::abs(Emin_), std::abs(Emax_)});
    trivial_case_ = std::abs(deltaE_) <= EPS * scale_ref;

    if (!trivial_case_) {
        norm_factor_ = 2.0 / deltaE_;
        x_ = (deltaE_ * cfg_.dt) / (2.0 * cfg_.hbar);
        phase_ = std::exp(I * ((0.5 * deltaE_ + Emin_) * cfg_.dt) / cfg_.hbar);
    } else {
        norm_factor_ = 0.0;
        x_ = 0.0;
        phase_ = {1.0, 0.0};
        trivial_phase_ = std::exp(I * (Emin_ * cfg_.dt) / cfg_.hbar);
    }

    const Eigen::Index n = T.a.size();
    p_prev_.resize(n);
    p_curr_.resize(n);
    tmp_.resize(n);
}

void ChebyshevEvolver::apply_normalized(const Eigen::VectorXcd& in,
                                        Eigen::VectorXcd& out) {
    tridiag_mul(T_, in, out);
    out.noalias() -= center_ * in;
    out *= norm_factor_;
}

StepResult ChebyshevEvolver::step(Eigen::VectorXcd& psi) {
    StepResult result;

    if (trivial_case_) {
        psi *= trivial_phase_;
        result.matvecs = 0;
        result.K_used = 1;
        result.bn_ratio = 1.0;
        return result;
    }

    if (p_prev_.size() != psi.size()) {
        const Eigen::Index n = psi.size();
        p_prev_.resize(n);
        p_curr_.resize(n);
        tmp_.resize(n);
    }

    const int nmax = std::max(0, cfg_.K);

    p_prev_ = psi;

    const double J0 = std::cyl_bessel_j(0, x_);
    Eigen::VectorXcd accum = J0 * p_prev_;

    double bn_abs_max = std::abs(J0);
    if (bn_abs_max == 0.0) {
        bn_abs_max = 0.0;
    }

    int terms = 1;
    double last_ratio = 1.0;
    std::complex<double> i_pow(1.0, 0.0);

    for (int n = 1; n <= nmax; ++n) {
        if (n == 1) {
            apply_normalized(p_prev_, p_curr_);
        } else {
            apply_normalized(p_curr_, tmp_);
            tmp_ *= 2.0;
            tmp_.noalias() -= p_prev_;
            p_prev_.swap(p_curr_);
            p_curr_.swap(tmp_);
        }

        i_pow *= I;
        const double Jn = std::cyl_bessel_j(n, x_);
        const std::complex<double> bn = 2.0 * i_pow * Jn;
        const double abs_bn = std::abs(bn);

        double ratio = 1.0;
        if (bn_abs_max > 0.0) {
            ratio = abs_bn / bn_abs_max;
        }

        accum.noalias() += bn * p_curr_;
        ++terms;

        if (abs_bn > bn_abs_max) {
            bn_abs_max = abs_bn;
        }

        last_ratio = (bn_abs_max > 0.0) ? ratio : 0.0;

        if (cfg_.tolerance > 0.0 && bn_abs_max > 0.0 && ratio < cfg_.tolerance) {
            break;
        }
    }

    psi = phase_ * accum;

    result.matvecs = std::max(0, terms - 1);
    result.K_used = terms;
    result.bn_ratio = last_ratio;
    return result;
}

