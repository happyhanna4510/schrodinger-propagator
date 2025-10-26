#pragma once

#include <Eigen/Core>

#include <iosfwd>
#include <optional>
#include <string>

#include "core/spectral.hpp"

struct StepMetrics {
    double norm = 0.0;
    double norm_err = 0.0;
    double theta = 0.0;
};

StepMetrics compute_step_metrics(const SpectralData& spectral,
                                 const Eigen::VectorXcd& psi,
                                 double t,
                                 double dx);

void write_step_csv_header(std::ofstream& f, bool include_cheb_extras);

void write_step_csv_row(std::ofstream& f,
                        const std::string& method,
                        int step,
                        double t,
                        double dt,
                        double dt_ms,
                        int matvecs,
                        const StepMetrics& metrics,
                        std::optional<int> K_used,
                        std::optional<double> bn_ratio,
                        bool include_cheb_extras);

void print_step_console(const std::string& method,
                        int step,
                        double t,
                        double dt_ms,
                        int matvecs,
                        const StepMetrics& metrics,
                        std::optional<int> K_used);

