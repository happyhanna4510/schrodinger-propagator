#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <sstream>
#include "grid.hpp"
#include "cli.hpp"
#include "snapshots.hpp"

#include <Eigen/Core>

namespace fs = std::filesystem;

// Печать таблицы энергий: num vs analytic
void print_energy_table(const std::vector<double>& Eval,
                        const std::vector<double>& Eanal,
                        int first);

// Проекция psi на первые k собственных векторов
std::vector<double> populations(const Eigen::MatrixXd& EV,
                                const Eigen::VectorXd& psi,
                                int k, double dx);

// --- CSV writers: принимают fs::path ---
void save_xy_csv(const fs::path& path,
                 const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::string& colx="x",
                 const std::string& coly="y");

void save_matrix_csv(const fs::path& path,
                     const std::vector<double>& x,
                     const Eigen::MatrixXd& cols,
                     int first_cols,
                     const std::string& colx="x",
                     const std::string& base="psi_n");

void save_vector_csv(const fs::path& path,
                     const std::vector<double>& v,
                     const std::string& col="value");



namespace io {
namespace fs = std::filesystem;

// ТОЛЬКО объявление (без inline и без тела)
fs::path make_csv_path(const fs::path& out_dir,
                       const Params& P,
                       int K, double dt,
                       const Grid& g);

// Объявление класса (без тел методов)
class WideDump {
public:
    WideDump(const fs::path& csv_path,
             const std::vector<double>* x_inner,
             bool write_re = false,
             bool write_im = false);

    bool enabled() const { return x_inner_ != nullptr; }
    void write(const Eigen::VectorXcd& psi, double t);

private:
    const std::vector<double>* x_inner_{nullptr};
    bool write_re_{false}, write_im_{false};

    fs::path parent_;
    std::string stem_;
    fs::path p_abs_, p_re_, p_im_;
    std::ofstream f_abs_;
    std::ofstream f_re_;
    std::ofstream f_im_;
    bool wrote_abs_ = false;
    bool wrote_re_ = false;
    bool wrote_im_ = false;
};
} // namespace io
