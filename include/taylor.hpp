#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex>
#include <filesystem>
#include <string>
#include <vector>

#include "cli.hpp"
#include "grid.hpp"

using cplx = std::complex<double>;

// ===== трёхдиагональная форма H =====
struct Tridiag {
    Eigen::VectorXd a; // main   (size M)
    Eigen::VectorXd b; // upper  (size M-1)
    Eigen::VectorXd c; // lower  (size M-1)
};

Tridiag make_tridiag_from_dense(const Eigen::Ref<const Eigen::MatrixXd>& H);

enum class LogExtras : unsigned {
    None   = 0u,
    P0     = 1u << 0,
    ErrPhi = 1u << 1,
    Both   = static_cast<unsigned>(P0) | static_cast<unsigned>(ErrPhi)
};

constexpr inline bool has(LogExtras e, LogExtras bit) {
    return (static_cast<unsigned>(e) & static_cast<unsigned>(bit)) != 0u;
}

// Запускает эволюцию Тейлора и пишет логи/wide.
// Возвращает путь к файлу лога CSV.
std::filesystem::path run_taylor_evolution(const Grid& g,
                                           const std::vector<double>& U_true,
                                           const Params& P,
                                           const std::filesystem::path& out_dir);

// Возвращаем простые метрики теста обратимости
struct RevTestResult {
    double err_back;     // ||psi_back - psi_init||2
    double norm0;        // норма начального состояния (sqrt)
    double norm_fwd;     // норма после прямого прогона (sqrt)
    double norm_back;    // норма после обратного прогона (sqrt)
    double max_drift;    // максимум |norm-1|, наблюдавшийся в обоих проходах
};

RevTestResult test_time_reversibility(const Tridiag& T,
                                      const Eigen::VectorXcd& psi_init,
                                      double dx, double dt, int nsteps, int K);
