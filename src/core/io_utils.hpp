#pragma once

#include <Eigen/Core>

#include <iosfwd>

#include "taylor.hpp"

#include "spectral.hpp"

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
