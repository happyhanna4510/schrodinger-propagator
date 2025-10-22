#include "snapshots.hpp"
#include <fstream>
#include <iomanip>

void save_psi_csv(const std::filesystem::path& path,
                  const std::vector<double>& x_inner,
                  const Eigen::VectorXcd& psi,
                  double t)
{
    std::ofstream f(path.string(), std::ios::out);
    f << std::setprecision(17);
    f << "t," << t << "\n";
    f << "x,RePsi,ImPsi,Abs2\n";
    const int M = (int)psi.size();
    for (int i=0; i<M; ++i){
        double re = psi[i].real();
        double im = psi[i].imag();
        double ab2 = re*re + im*im;
        f << x_inner[i] << "," << re << "," << im << "," << ab2 << "\n";
    }
}


void append_abs2_wide(const std::filesystem::path& path,
                      const std::vector<double>& x,
                      const Eigen::VectorXcd& psi,
                      double t,
                      bool write_header)
{
    std::ofstream f(path.string(), std::ios::out | std::ios::app);
    f << std::setprecision(16);

    if (write_header) {
        f << "t";
        for (double xi : x) f << "," << xi;
        f << "\n";
    }

    f << t;
    for (int i = 0; i < psi.size(); ++i) {
        double re = psi[i].real(), im = psi[i].imag();
        double abs2 = re*re + im*im;
        f << "," << abs2;        // хочешь — замени на re/im
    }
    f << "\n";
}

void append_re_wide (const std::filesystem::path& path,
                     const std::vector<double>& x,
                     const Eigen::VectorXcd& psi, double t,
                     bool write_header){
    std::ofstream f(path.string(), std::ios::out | std::ios::app);
    f << std::setprecision(16);
    if (write_header){ f << "t"; for (double xi: x) f << "," << xi; f << "\n"; }
    f << t;
    for (int i=0;i<psi.size();++i) f << "," << psi[i].real();
    f << "\n";
}

void append_im_wide (const std::filesystem::path& path,
                     const std::vector<double>& x,
                     const Eigen::VectorXcd& psi, double t,
                     bool write_header){
    std::ofstream f(path.string(), std::ios::out | std::ios::app);
    f << std::setprecision(16);
    if (write_header){ f << "t"; for (double xi: x) f << "," << xi; f << "\n"; }
    f << t;
    for (int i=0;i<psi.size();++i) f << "," << psi[i].imag();
    f << "\n";
}
