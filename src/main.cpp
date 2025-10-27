#include <filesystem>
#include <iostream>
#if defined(_MSC_VER) || defined(__SSE__)
#    include <xmmintrin.h>
#endif
#if defined(_MSC_VER) || defined(__SSE3__)
#    include <pmmintrin.h>
#endif

#include "cli.hpp"
#include "grid.hpp"
#include "morse_potential.hpp"
#include "morse_static.hpp"
#include "paths.hpp"
#include "runtime_evolution.hpp"

namespace fs = std::filesystem;

int main(int argc, char** argv) 
{

#if defined(_MSC_VER) || defined(__SSE__)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
#if defined(_MSC_VER) || defined(__SSE3__)
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    Params P = parse_args(argc, argv);

    Grid g(P.N, P.xmax);
    auto U_true = morse_potential(g, P.gamma);

    fs::path out = resolve_outdir(fs::path(argv[0]), P.outdir);
    fs::create_directories(out);

    if (!P.evolve_only) {
        write_morse_outputs(g, U_true, P.gamma, P.first, out, P.quiet);
    } else if (!P.quiet) {
        std::cout << "# --evolve_only: skipping Morse static outputs.\n";
    }

    if (!P.do_evolve) {
        if (!P.quiet) {
            std::cout << "# No evolution (use --evolve [taylor|rk4])\n";
        }
        return 0;
    }

    
    fs::path csv_path = run_time_evolution(g, U_true, P, out);

    return 0;
}

