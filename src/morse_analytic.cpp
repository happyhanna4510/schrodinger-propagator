#include "morse_analytic.hpp"
#include <vector>
#include <cmath>
#include <algorithm>

int morse_nmax(double gamma){
    // максимально допустимый n для связанных уровней
    int nmax = (int)std::floor(gamma - 0.5 - 1e-12);
    return std::max(nmax, -1); // -1 если gamma<=0.5 (нет связанных)
}

std::vector<double> morse_energies(double gamma, int how_many){
    int nmax = morse_nmax(gamma);
    int m = std::min(how_many, std::max(nmax+1, 0));
    std::vector<double> E; E.reserve(m);
    for (int n=0;n<m;n++){
        double En = 1.0 - std::pow(gamma - (n + 0.5), 2) / (gamma*gamma);
        E.push_back(En);
    }
    return E;
}
