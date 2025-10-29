#include "paths.hpp"

fs::path find_project_root(const fs::path& exe_dir) 
{
    fs::path cur = exe_dir;
    for (int i = 0; i < 10; ++i) {
        if (exists(cur / ".git") || exists(cur / "CMakeLists.txt")) return cur;
        if (!cur.has_parent_path()) break;
        cur = cur.parent_path();
    }
    return exe_dir;
}
