#pragma once
#include "grid.hpp"
#include <vector>
#include <filesystem>

void write_morse_outputs(const Grid& g,
                                const std::vector<double>& U_true,
                                int gamma, int first,
                                const std::filesystem::path& out);
