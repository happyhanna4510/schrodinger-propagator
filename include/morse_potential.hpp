#pragma once
#include <vector>
#include "grid.hpp"

// Морзе в безразмерной форме с параметром gamma:
// U(x) = 1 - ( (gamma - n - 1/2)^2 / gamma^2 ) даёт аналитические уровни,
// а сам профиль U(x) берём как стандартный: U(x) = (1 - e^{-x})^2
// (gamma влияет на спектр — см. morse_analytic.hpp)

std::vector<double> morse_potential(const Grid& g, double gamma);

// (опционально) жёсткое отсечение сверху для численной стабильности:
inline void cap_potential(std::vector<double>& U, double Vmax) {
    for (auto& v : U) if (v > Vmax) v = Vmax;
}