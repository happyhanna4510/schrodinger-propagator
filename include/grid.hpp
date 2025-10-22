#pragma once
#include <vector>

struct Grid {
    int    N;       // число узлов, включая края
    double xmax;    // правая граница, левая = -xmax
    double dx;      // шаг
    std::vector<double> x; // координаты

    Grid(int N_, double xmax_);
};
