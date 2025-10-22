#include "hamiltonian.hpp"
#include <cassert>

Eigen::MatrixXd build_hamiltonian(const Grid& g, const std::vector<double>& U){
    assert((int)U.size()==g.N);
    const int M = g.N - 2;          // внутренние точки
    const double dx = g.dx;
    const double invdx2 = 1.0/(dx*dx);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, M);
    // Кинетика: -1/2 d^2/dx^2 ≈ diag(+1/dx^2), offdiag(-1/(2 dx^2))
    const double diagK  = 1.0 * invdx2;
    const double offK   = -0.5 * invdx2;

    for (int i=0;i<M;i++){
        H(i,i) = diagK + U[i+1];     // + U на соответствующей внутренней точке
        if (i>0)   H(i, i-1) = offK;
        if (i<M-1) H(i, i+1) = offK;
    }
    return H;
}

void renormalize_eigenvectors(Eigen::MatrixXd& EV, double dx){
    const int M = (int)EV.rows();
    const int K = (int)EV.cols();
    for (int k=0;k<K;k++){
        double norm = 0.0;
        for (int i=0;i<M;i++) norm += EV(i,k)*EV(i,k);
        norm *= dx;
        EV.col(k) /= std::sqrt(norm);
    }
}
