#include "io.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <sstream>

using std::cout; using std::setw; using std::setprecision; using std::fixed;

void print_energy_table(const std::vector<double>& Eval,
                        const std::vector<double>& Eanal, int first) {
    cout << " n  " << setw(16) << "E_num" << setw(16) << "E_anal"
         << setw(16) << "abs diff" << "\n";
    int m = std::min((int)Eval.size(), (int)Eanal.size());
    for (int n=0; n<std::min(first,m); ++n){
        double d = std::abs(Eval[n]-Eanal[n]);
        cout << setw(2) << n
             << setw(16) << setprecision(10) << fixed << Eval[n]
             << setw(16) << setprecision(10) << fixed << Eanal[n]
             << setw(16) << std::scientific << setprecision(3) << d << "\n";
    }
}

std::vector<double> populations(const Eigen::MatrixXd& EV,
                                const Eigen::VectorXd& psi, int k, double dx) {
    int K = std::min(k, (int)EV.cols());
    std::vector<double> pop(K,0.0);
    for (int j=0;j<K;j++){
        double cj = (EV.col(j).array() * psi.array()).sum() * dx;
        pop[j] = cj*cj; // все реально => |c|^2 = c^2
    }
    return pop;
}

// --- CSV writers on wide-safe paths ---
static std::ofstream open_ofs(const fs::path& p) {
    std::ofstream f(p, std::ios::binary);
    if (!f) {
        std::wcerr << L"[io] cannot open for write: " << p.wstring() << std::endl;
    }
    return f;
}

void save_xy_csv(const fs::path& path,
                 const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::string& colx,
                 const std::string& coly) {
    std::ofstream f = open_ofs(path);
    if (!f) return;
    f << colx << "," << coly << "\n";
    size_t n = std::min(x.size(), y.size());
    for (size_t i=0;i<n;i++) f << std::setprecision(12) << x[i] << "," << y[i] << "\n";
}

void save_matrix_csv(const fs::path& path,
                     const std::vector<double>& x,
                     const Eigen::MatrixXd& cols,
                     int first_cols,
                     const std::string& colx,
                     const std::string& base) {
    int M = (int)x.size();
    int K = std::min(first_cols, (int)cols.cols());
    std::ofstream f = open_ofs(path);
    if (!f) return;
    f << colx;
    for (int j=0;j<K;j++) f << "," << base << j;
    f << "\n";
    for (int i=0;i<M;i++){
        f << std::setprecision(12) << x[i];
        for (int j=0;j<K;j++) f << "," << std::setprecision(12) << cols(i,j);
        f << "\n";
    }
}

void save_vector_csv(const fs::path& path,
                     const std::vector<double>& v,
                     const std::string& col) {
    std::ofstream f = open_ofs(path);
    if (!f) return;
    f << col << "\n";
    for (auto& val : v) f << std::setprecision(12) << val << "\n";
}




static std::string fmt_dt_for_filename(double dt) {
    // считаем, что dt > 0; при нуле вернём "0"
    if (dt == 0.0) return "0";

    // проверка "dt — степень десяти" с малой погрешностью двойной точности
    const double log10v = std::log10(dt);
    const long long e   = std::llround(log10v);
    const double pow10e = std::pow(10.0, (double)e);
    const double relerr = std::fabs(dt - pow10e) / dt;

    if (relerr < 1e-12) {
        // это точно 10^e → печатаем "1e-4" (без + и ведущих нулей)
        if (e == 0) return "1";
        return std::string("1e") + (e < 0 ? "" : "") + std::to_string(e);
        // для положительных экспонент получится "1e3", для отрицательных — "1e-4"
    }

    // иначе — компактно, как у тебя
    std::ostringstream os;
    os << std::setprecision(6) << std::defaultfloat << dt;
    return os.str();
}


namespace io {
namespace fs = std::filesystem;

// РЕАЛИЗАЦИЯ (без inline)
fs::path make_csv_path(const fs::path& out_dir,
                       const Params& P,
                       int K, double dt,
                       const Grid& g)
{
    fs::create_directories(out_dir);

    auto fmtg = [](double v){
        std::ostringstream os;
        os << std::setprecision(6) << std::defaultfloat << v; // компактно
        return os.str();
    };

    return !P.csv_name.empty()
        ? out_dir / P.csv_name
        : out_dir / ("taylor_K" + std::to_string(K) +
                     "_dt" + fmt_dt_for_filename(dt) +
                     "_g"  + fmtg(P.gamma) +
                     "_N"  + std::to_string(g.N) +
                     "_x"  + fmtg(g.xmax) + ".csv");
}

// Конструктор WideDump — РЕАЛИЗАЦИЯ
WideDump::WideDump(const fs::path& csv_path,
                   const std::vector<double>* x_inner,
                   bool write_re, bool write_im)
: x_inner_(x_inner), write_re_(write_re), write_im_(write_im)
{
    if (!x_inner_) return;
    parent_ = csv_path.parent_path();
    stem_   = csv_path.stem().string();
    p_abs_  = parent_ / (stem_ + "_abs2_wide.csv");
    if (write_re_) p_re_ = parent_ / (stem_ + "_re_wide.csv");
    if (write_im_) p_im_ = parent_ / (stem_ + "_im_wide.csv");
}

// Метод write — РЕАЛИЗАЦИЯ
void WideDump::write(const Eigen::VectorXcd& psi, double t) {
    if (!enabled()) return;

    ::append_abs2_wide(p_abs_, *x_inner_, psi, t, !wrote_abs_);
    wrote_abs_ = true;

    if (write_re_) {
        ::append_re_wide(p_re_, *x_inner_, psi, t, !wrote_re_);
        wrote_re_ = true;
    }
    if (write_im_) {
        ::append_im_wide(p_im_, *x_inner_, psi, t, !wrote_im_);
        wrote_im_ = true;
    }
}

} // namespace io
