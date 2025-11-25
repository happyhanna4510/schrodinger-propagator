#include "evolve/evolve_chebyshev.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include "core/math_utils.hpp"

namespace 
{ //i = sqrt(-1)
constexpr std::complex<double> I(0.0, 1.0);
constexpr double EPS = 1e-14;
}

namespace fs = std::filesystem;

ChebyshevEvolver::ChebyshevEvolver(const Tridiag& T,
                                   const SpectralData& spectral,
                                   const EvolverConfig& cfg)
    : EvolverBase(T, cfg) 
{
    //ekstrema widma min i maks energia
    Emin_ = spectral.evals.minCoeff();
    Emax_ = spectral.evals.maxCoeff();
    //szerokość przedziału widma i jego środek
    deltaE_ = Emax_ - Emin_;
    center_ = 0.5 * (Emax_ + Emin_);

    //gdy widmo jest praktycznie punktowe
    const double scale_ref = std::max({1.0, std::abs(Emin_), std::abs(Emax_)});
    trivial_case_ = std::abs(deltaE_) <= EPS * scale_ref;

    if (!trivial_case_) 
    {//współczynnik normalizacji do zakresu [-1, 1] dla wielomianów Czebyszewa
        norm_factor_ = 2.0 / deltaE_;
        // x dla współczynników Bessela J_n(x) 
        x_ = (deltaE_ * cfg_.dt) / (2.0 * cfg_.hbar);
        //globalna faza
        phase_ = std::exp(-I * ((0.5 * deltaE_ + Emin_) * cfg_.dt) / cfg_.hbar);

    } else 
    {
        norm_factor_ = 0.0;
        x_ = 0.0;
        phase_ = {1.0, 0.0};
        trivial_phase_ = std::exp(-I * (Emin_ * cfg_.dt) / cfg_.hbar);

    }

    const Eigen::Index n = T.a.size();
    p_prev_.resize(n);
    p_curr_.resize(n);
    tmp_.resize(n);

    cheb_beta_log_path_ = cfg.cheb_beta_log_path;
    log_betas_ = !cheb_beta_log_path_.empty();
}

void ChebyshevEvolver::apply_normalized(const Eigen::VectorXcd& in,
                                        Eigen::VectorXcd& out)
{
    tridiag_mul(T_, in, out);
    out.noalias() -= center_ * in;
    out *= norm_factor_;
}

void ChebyshevEvolver::maybe_write_betas(int terms)
{
    if (!log_betas_ || betas_written_) {
        return;
    }

    const int available = static_cast<int>(betas_.size());
    const int count = std::min(terms, available);

    fs::path parent = cheb_beta_log_path_.parent_path();
    std::error_code ec;
    if (!parent.empty()) {
        fs::create_directories(parent, ec);
    }

    std::ofstream out(cheb_beta_log_path_, std::ios::out | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("failed to open Chebyshev beta log file");
    }
    out << std::setprecision(16);
    out << "n,beta_real,beta_imag,beta_abs\n";

    for (int n = 0; n < count; ++n) {
        const std::complex<double>& beta = betas_[static_cast<std::size_t>(n)];
        out << n << ',' << beta.real() << ',' << beta.imag() << ','
            << std::abs(beta) << '\n';
    }

    betas_written_ = true;
}

StepResult ChebyshevEvolver::step(Eigen::VectorXcd& psi)
{
    StepResult result;

    betas_.clear();

    if (trivial_case_) //bez matveców
    {
        psi *= trivial_phase_;
        const double J0 = std::cyl_bessel_j(0, x_);
        betas_.push_back(std::complex<double>(J0, 0.0));
        maybe_write_betas(1);
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

    const int max_safe = 10000; //zabezpieczenie przed zapętleniem

    p_prev_ = psi;

    betas_.reserve(static_cast<std::size_t>(std::max(2, cfg_.K)));

    const double J0 = std::cyl_bessel_j(0, x_);
    betas_.push_back(std::complex<double>(J0, 0.0));
    Eigen::VectorXcd accum = J0 * p_prev_;

    double bn_abs_max = std::abs(J0);
    int terms = 1;
    double last_ratio = 1.0;
    std::complex<double> i_pow(1.0, 0.0); // (-i)^0

    for (int n = 1; n < max_safe; ++n) {
        //rekurencja Czebyszewa
        if (n == 1) {
            apply_normalized(p_prev_, p_curr_);
        } else {
            apply_normalized(p_curr_, tmp_);
            tmp_ *= 2.0;
            tmp_.noalias() -= p_prev_;
            p_prev_.swap(p_curr_);
            p_curr_.swap(tmp_);
        }

        i_pow *= -I; // (-i)^n

        const double Jn = std::cyl_bessel_j(n, x_);
        const std::complex<double> bn = 2.0 * i_pow * Jn;
        const double abs_bn = std::abs(bn);

        betas_.push_back(bn);

        if (abs_bn > bn_abs_max) {
            bn_abs_max = abs_bn;
        }

        double ratio = (bn_abs_max > 0.0) ? abs_bn / bn_abs_max : 0.0;
        last_ratio = ratio;

        accum.noalias() += bn * p_curr_;
        ++terms;

        //ratio < tolerance -przerwij
        if (cfg_.tolerance > 0.0 && ratio < cfg_.tolerance) {
            break;
        }
    }


    psi = phase_ * accum;
    maybe_write_betas(terms);
    //liczba matveców
    result.matvecs = std::max(0, terms - 1);
    result.K_used = terms;
    result.bn_ratio = last_ratio;
    return result;
}

