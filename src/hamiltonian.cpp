#include "hamiltonian.hpp"
#include <cassert>

Eigen::MatrixXd build_hamiltonian(const Grid& g, const std::vector<double>& U)
{
    assert((int)U.size()==g.N);
    const int M = g.N - 2;          // liczba punktów wewnętrznych (bez brzegów)
    const double dx = g.dx;
    const double invdx2 = 1.0/(dx*dx);

    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, M);
    //  (schemat 3-punktowy)
    const double diagK  = 1.0 * invdx2;     // wkład diagonalny z laplasjanu
    const double offK   = -0.5 * invdx2;    // sąsiedzi (i-1, i+1)

    for (int i=0;i<M;i++)
    {
        H(i,i) = diagK + U[i+1];     // + U  odpowiadającym punkcie wewnętrznym
        if (i>0)   H(i, i-1) = offK;
        if (i<M-1) H(i, i+1) = offK;
    }
    return H;
}

//normalizacja L2 na siatce
void renormalize_eigenvectors(Eigen::MatrixXd& EV, double dx)
{
    const int M = (int)EV.rows();
    const int K = (int)EV.cols();
    for (int k=0;k<K;k++)
    {
        double norm = 0.0;
        for (int i=0;i<M;i++) norm += EV(i,k)*EV(i,k);
        norm *= dx;
        EV.col(k) /= std::sqrt(norm);
    }
}
