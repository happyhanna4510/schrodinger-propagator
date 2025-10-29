#pragma once
#include <vector>
#include "grid.hpp"

//Potencjał Morse’a na siatce: U(x) = (1 - e^{-a x})^2,

std::vector<double> morse_potential(const Grid& g, double gamma);

// ograniczenia na potencjał: U(x) = min(U(x), Vmax)
inline void cap_potential(std::vector<double>& U, double Vmax) {
    for (auto& v : U) if (v > Vmax) v = Vmax;
}