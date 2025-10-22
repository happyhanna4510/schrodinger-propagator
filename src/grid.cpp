#include "grid.hpp"

Grid::Grid(int N_, double xmax_) : N(N_), xmax(xmax_), dx(2.0*xmax_/ (N_-1)), x(N_) {
    for (int i=0;i<N;i++) x[i] = -xmax + dx*i;
}
