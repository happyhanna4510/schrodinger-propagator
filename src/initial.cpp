#include "initial.hpp"
#include <cmath>

Eigen::VectorXd gaussian_on_inner(const Grid& g, double x0, double sigma){
    const int M = g.N - 2;
    Eigen::VectorXd v(M);
    for (int i=0;i<M;i++){
        double x = g.x[i+1];
        v(i) = std::exp(-(x-x0)*(x-x0)/(2*sigma*sigma));
    }
    double norm = (v.array().square().sum()) * g.dx;
    v /= std::sqrt(norm);
    return v;
}

// complex Gauss with plane-wave phase e^{i k0 x} (удобно для "движущейся" пачки)
using cplx = std::complex<double>;
Eigen::VectorXcd gaussian_complex_on_inner(const Grid& g, double x0, double sigma, double k0){
    const int M = g.N - 2;
    Eigen::VectorXcd v(M);
    for (int i=0;i<M;i++){
        double x = g.x[i+1];
        double amp = std::exp(-(x-x0)*(x-x0)/(2*sigma*sigma));
        v(i) = cplx(amp, 0.0) * std::exp(cplx(0.0, k0 * x));
    }
    double S = 0.0; for(int i=0;i<M;++i) S += std::norm(v(i));
    v /= std::sqrt(S * g.dx);
    return v;
}
