#include "snapshots.hpp"

#include <fstream>
#include <iomanip>

void save_psi_csv(const std::filesystem::path& path,
                  const std::vector<double>& x_inner,
                  const Eigen::VectorXcd& psi,
                  double t) {
    std::ofstream f(path, std::ios::out | std::ios::trunc);
    if (!f) {
        return;
    }

    f << std::setprecision(17);
    f << "t," << t << "\n";
    f << "x,RePsi,ImPsi,Abs2\n";

    const int M = static_cast<int>(psi.size());
    for (int i = 0; i < M; ++i) {
        const double re = psi[i].real();
        const double im = psi[i].imag();
        const double ab2 = re * re + im * im;
        f << x_inner[i] << ',' << re << ',' << im << ',' << ab2 << "\n";
    }
}

void append_abs2_wide(const std::filesystem::path& path,
                      const std::vector<double>& x,
                      const Eigen::VectorXcd& psi,
                      double t,
                      bool write_header) {
    std::ofstream f(path, std::ios::out | std::ios::app);
    if (!f) {
        return;
    }

    f << std::setprecision(16);
    if (write_header) {
        f << "t";
        for (double xi : x) f << ',' << xi;
        f << "\n";
    }

    f << t;
    for (int i = 0; i < psi.size(); ++i) {
        const double re = psi[i].real();
        const double im = psi[i].imag();
        f << ',' << (re * re + im * im);
    }
    f << "\n";
}

void append_re_wide(const std::filesystem::path& path,
                    const std::vector<double>& x,
                    const Eigen::VectorXcd& psi,
                    double t,
                    bool write_header) {
    std::ofstream f(path, std::ios::out | std::ios::app);
    if (!f) {
        return;
    }

    f << std::setprecision(16);
    if (write_header) {
        f << "t";
        for (double xi : x) f << ',' << xi;
        f << "\n";
    }

    f << t;
    for (int i = 0; i < psi.size(); ++i) {
        f << ',' << psi[i].real();
    }
    f << "\n";
}

void append_im_wide(const std::filesystem::path& path,
                    const std::vector<double>& x,
                    const Eigen::VectorXcd& psi,
                    double t,
                    bool write_header) {
    std::ofstream f(path, std::ios::out | std::ios::app);
    if (!f) {
        return;
    }

    f << std::setprecision(16);
    if (write_header) {
        f << "t";
        for (double xi : x) f << ',' << xi;
        f << "\n";
    }

    f << t;
    for (int i = 0; i < psi.size(); ++i) {
        f << ',' << psi[i].imag();
    }
    f << "\n";
}

