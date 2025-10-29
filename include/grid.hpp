#pragma once
#include <vector>

struct Grid {
    int    N;       // liczba punktów siatki
    double xmax;    // granica x w obu kierunkach
    double dx;      // krok siatki
    std::vector<double> x; // koordynaty punktów siatki

    Grid(int N_, double xmax_);
};
