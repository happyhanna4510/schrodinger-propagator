#pragma once
#include <Eigen/Dense>
#include <filesystem>
#include <vector>

void save_psi_csv(const std::filesystem::path& path,
                  const std::vector<double>& x_inner,
                  const Eigen::VectorXcd& psi,
                  double t);


void append_abs2_wide(const std::filesystem::path& path,
                      const std::vector<double>& x,
                      const Eigen::VectorXcd& psi,
                      double t,
                      bool write_header);

                      void append_im_wide (const std::filesystem::path& path,
                     const std::vector<double>& x,
                     const Eigen::VectorXcd& psi, double t,
                     bool write_header);

void append_re_wide (const std::filesystem::path& path,
                     const std::vector<double>& x,
                     const Eigen::VectorXcd& psi, double t,
                     bool write_header);