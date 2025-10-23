#pragma once

#include <filesystem>
#include <vector>

#include "cli.hpp"
#include "grid.hpp"

namespace fs = std::filesystem;

fs::path run_time_evolution(const Grid& g,
                            const std::vector<double>& U_true,
                            const Params& P,
                            const fs::path& out_dir);

fs::path run_taylor_evolution(const Grid& g,
                              const std::vector<double>& U_true,
                              const Params& P,
                              const fs::path& out_dir);

