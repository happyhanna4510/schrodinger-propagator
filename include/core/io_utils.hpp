#pragma once

#include <Eigen/Core>

#include <iosfwd>

enum class LogExtras : unsigned {
    None   = 0u,
    P0     = 1u << 0,
    ErrPhi = 1u << 1,
    Both   = static_cast<unsigned>(P0) | static_cast<unsigned>(ErrPhi)
};

constexpr inline bool has(LogExtras e, LogExtras bit) {
    return (static_cast<unsigned>(e) & static_cast<unsigned>(bit)) != 0u;
}

struct SpectralData;

struct LogSnapshot {
    int    step = 0;
    double t = 0.0;
    double norm_sq = 0.0;
    double norm = 0.0;
    double drift_sq = 0.0;
    double drift = 0.0;
    double p0 = 0.0;
    double err_exact = 0.0;
    double err_true = 0.0;
    double theta = 0.0;
    double edge_left = 0.0;
    double edge_right = 0.0;
};

LogSnapshot collect_snapshot(const SpectralData& spectral,
                             const Eigen::VectorXcd& psi,
                             double t,
                             double dx,
                             int step,
                             double norm_sq0,
                             double norm0);

void write_csv_header(std::ofstream& f, LogExtras extras);
void write_csv_row(std::ofstream& f, const LogSnapshot& snap, LogExtras extras);
void print_snapshot(const LogSnapshot& snap, LogExtras extras);

