#pragma once
#include <Eigen/Dense>
#include <complex>
#include <string>
#include <filesystem>
#include <vector>
#include "grid.hpp"
#include "cli.hpp"
#include <Eigen/Core>
#include <complex>
#include <string>
#include <vector>


using cplx = std::complex<double>;

struct TaylorLogRow {
    int    step;
    double t;
    double norm;
    double p0;
    double err_exact;
};

// ===== трёхдиагональная форма H =====
struct Tridiag {
    Eigen::VectorXd a; // main   (size M)
    Eigen::VectorXd b; // upper  (size M-1)
    Eigen::VectorXd c; // lower  (size M-1)
};

Tridiag make_tridiag_from_dense(const Eigen::MatrixXd& H);

// ===== метрики/нормы (без изменений) =====
double normL2_dx(const Eigen::VectorXcd& psi, double dx);

TaylorLogRow metrics_vs_ground(const Eigen::VectorXcd& psi,
                               const Eigen::VectorXcd& phi0,
                               double E0, double t, double dx,
                               int step);

// ====== быстрый тейлор через трёхдиаг ======
void taylor_step_tridiag(const Tridiag& T, Eigen::VectorXcd& psi, double dt, int K);


enum class LogExtras : unsigned {
    None   = 0u,
    P0     = 1u << 0,
    ErrPhi = 1u << 1,
    Both   = static_cast<unsigned>(P0) | static_cast<unsigned>(ErrPhi) // <- cast
};

constexpr inline bool has(LogExtras e, LogExtras bit) {
    return (static_cast<unsigned>(e) & static_cast<unsigned>(bit)) != 0u;
}

// Объявление функции (вариант с параметром по умолчанию)
void evolve_taylor_tridiag(const Tridiag& T,
                           const Eigen::VectorXcd& psi_init,
                           const Eigen::VectorXcd& phi0,
                           double E0, double dx, double dt, int nsteps, int K,
                           int log_every, const std::string& csv_path,
                           const std::vector<double>* x_inner,
                           const std::string& snaps_dir,
                           LogExtras extras = LogExtras::None);


// Запускает эволюцию Тейлора и пишет логи/wide.
// Возвращает путь к файлу лога CSV.
std::filesystem::path run_taylor_evolution(const Grid& g,
                                           const std::vector<double>& U_true,
                                           const Params& P,
                                           const std::filesystem::path& out_dir);


                                           
// Возвращаем простые метрики теста обратимости
struct RevTestResult {
    double err_back;     // ||psi_back - psi_init||2
    double norm0;        // норма начального состояния
    double norm_fwd;     // норма после прямого прогона
    double norm_back;    // норма после обратного прогона
    double max_drift;    // максимум |norm-1|, наблюдавшийся в обоих проходах
};

RevTestResult test_time_reversibility(const Tridiag& T,
                                      const Eigen::VectorXcd& psi_init,
                                      double dx, double dt, int nsteps, int K);
