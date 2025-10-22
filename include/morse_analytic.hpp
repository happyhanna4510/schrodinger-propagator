#pragma once
#include <vector>

// Аналитические энергии Морзе в принятой безразмерной форме:
// E_n = 1 - ((gamma - n - 1/2)^2)/(gamma^2),   n = 0,1,..., n_max
// где n_max = floor(gamma - 1/2 - eps)
int morse_nmax(double gamma);
std::vector<double> morse_energies(double gamma, int how_many);
