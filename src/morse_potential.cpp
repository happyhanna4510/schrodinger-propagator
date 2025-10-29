#include "morse_potential.hpp"
#include <algorithm>
#include <cmath>

// Potencjał Morse’a na siatce: U(x) = (1 - e^{-a x})^2,
// gdzie a = sqrt(2)/gamma 
std::vector<double> morse_potential(const Grid& g, double gamma) 
{
    std::vector<double> U(g.N);
    const double a = std::sqrt(2.0) / gamma;   
    for (int i=0;i<g.N;i++) {
        double ex = std::exp(-a * g.x[i]);
        U[i] = (1.0 - ex) * (1.0 - ex);        
    }
    return U;
}

