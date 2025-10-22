#pragma once
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

fs::path find_project_root(const fs::path& exe_dir);

inline fs::path resolve_outdir(const fs::path& exe_path, const std::string& outdir){
    fs::path exe_dir = fs::absolute(exe_path).parent_path();
    fs::path root = find_project_root(exe_dir);
    fs::path out = fs::path(outdir);
    if (!out.is_absolute()) out = root / out;
    return out;
}
