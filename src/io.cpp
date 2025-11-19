#include "io.hpp"
#include <iostream>
#include <iomanip>

#include <fstream>
#include <cmath>

#include <filesystem>
#include <sstream>

using std::cout; using std::setw; using std::setprecision; using std::fixed;

void print_energy_table(const std::vector<double>& Eval,
                        const std::vector<double>& Eanal, int first) 
{
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
                                const Eigen::VectorXd& psi, int k, double dx) 
{
    int K = std::min(k, (int)EV.cols());
    std::vector<double> pop(K,0.0);
    for (int j=0;j<K;j++){
        double cj = (EV.col(j).array() * psi.array()).sum() * dx;
        pop[j] = cj*cj; // все реально => |c|^2 = c^2
    }
    return pop;
}

// --- CSV writers on wide-safe paths ---
static std::ofstream open_ofs(const fs::path& p) 
{
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
                 const std::string& coly) 
{
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
                     const std::string& base) 
{
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
                     const std::string& col) 
{
    std::ofstream f = open_ofs(path);
    if (!f) return;
    f << col << "\n";
    for (auto& val : v) f << std::setprecision(12) << val << "\n";
}




static std::string fmt_dt_for_filename(double dt) 
{
    if (dt == 0.0) return "0";

    const double log10v = std::log10(dt);
    const long long e   = std::llround(log10v);
    const double pow10e = std::pow(10.0, (double)e);
    const double relerr = std::fabs(dt - pow10e) / dt;

    if (relerr < 1e-12) {
        if (e == 0) return "1";
        return std::string("1e") + (e < 0 ? "" : "") + std::to_string(e);
    }

    std::ostringstream os;
    os << std::setprecision(6) << std::defaultfloat << dt;
    return os.str();
}


namespace io {
namespace fs = std::filesystem;

namespace {

void write_wide_header(std::ofstream& f, const std::vector<double>& x) {
    f << "t";
    for (double xi : x) {
        f << ',' << xi;
    }
    f << '\n';
}

void write_abs2_row(std::ofstream& f, const Eigen::VectorXcd& psi, double t) {
    f << t;
    for (int i = 0; i < psi.size(); ++i) {
        const double re = psi[i].real();
        const double im = psi[i].imag();
        f << ',' << (re * re + im * im);
    }
    f << '\n';
}

void write_component_row(std::ofstream& f, const Eigen::VectorXcd& psi, double t, bool real_part) {
    f << t;
    for (int i = 0; i < psi.size(); ++i) {
        f << ',' << (real_part ? psi[i].real() : psi[i].imag());
    }
    f << '\n';
}

} // namespace

fs::path make_csv_path(const fs::path& out_dir,
                       const Params& P,
                       const std::string& method,
                       int K, double dt,
                       const Grid& g)
{
    fs::create_directories(out_dir);

    auto fmtg = [](double v){
        std::ostringstream os;
        os << std::setprecision(6) << std::defaultfloat << v; 
        return os.str();
    };

    if (!P.csv_name.empty()) {
        return out_dir / P.csv_name;
    }

    std::ostringstream stem;
    stem << method;
    if (method == "taylor") {
        stem << "_K" << K;
    }
    stem << "_dt" << fmt_dt_for_filename(dt)
         << "_g"  << fmtg(P.gamma)
         << "_N"  << std::to_string(g.N)
         << "_x"  << fmtg(g.xmax);

    return out_dir / (stem.str() + ".csv");
}

fs::path make_energy_csv_path(const fs::path& log_csv_path)
{
    fs::path parent = log_csv_path.parent_path();
    if (parent.empty()) {
        parent = fs::current_path();
    }

    const std::string stem = log_csv_path.stem().string();
    fs::create_directories(parent);

    fs::path candidate = parent / (stem + "_energy.csv");
    int idx = 1;
    while (fs::exists(candidate)) {
        candidate = parent / (stem + "_energy_" + std::to_string(idx) + ".csv");
        ++idx;
    }
    return candidate;
}


WideDump::WideDump(const fs::path& csv_path,
                   const std::vector<double>* x_inner,
                   bool write_re, bool write_im)
: x_inner_(x_inner), write_re_(write_re), write_im_(write_im)
{
    if (!x_inner_) {
        return;
    }

    parent_ = csv_path.parent_path();
    if (!parent_.empty()) {
        std::error_code ec;
        fs::create_directories(parent_, ec);
    }

    stem_   = csv_path.stem().string();
    p_abs_  = parent_ / (stem_ + "_abs2_wide.csv");
    p_re_   = parent_ / (stem_ + "_re_wide.csv");
    p_im_   = parent_ / (stem_ + "_im_wide.csv");

    f_abs_.open(p_abs_, std::ios::out | std::ios::trunc);
    if (!f_abs_) {
        std::wcerr << L"[io] cannot open wide file: " << p_abs_.wstring() << std::endl;
        x_inner_ = nullptr;
        write_re_ = false;
        write_im_ = false;
        return;
    }
    f_abs_ << std::setprecision(16);

    if (write_re_) 
    {
        f_re_.open(p_re_, std::ios::out | std::ios::trunc);
        if (!f_re_) {
            std::wcerr << L"[io] cannot open wide file: " << p_re_.wstring() << std::endl;
            write_re_ = false;
        } else {
            f_re_ << std::setprecision(16);
        }
    }

    if (write_im_) 
    {
        f_im_.open(p_im_, std::ios::out | std::ios::trunc);
        if (!f_im_) {
            std::wcerr << L"[io] cannot open wide file: " << p_im_.wstring() << std::endl;
            write_im_ = false;
        } else {
            f_im_ << std::setprecision(16);
        }
    }
}

void WideDump::write(const Eigen::VectorXcd& psi, double t) 
{
    if (!enabled()) {
        return;
    }

    if (!wrote_abs_) {
        write_wide_header(f_abs_, *x_inner_);
        wrote_abs_ = true;
    }
    write_abs2_row(f_abs_, psi, t);

    if (write_re_) {
        if (!wrote_re_) {
            write_wide_header(f_re_, *x_inner_);
            wrote_re_ = true;
        }
        write_component_row(f_re_, psi, t, true);
    }

    if (write_im_) {
        if (!wrote_im_) {
            write_wide_header(f_im_, *x_inner_);
            wrote_im_ = true;
        }
        write_component_row(f_im_, psi, t, false);
    }
}

} // namespace io
