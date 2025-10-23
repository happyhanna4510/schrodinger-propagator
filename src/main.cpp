#include <filesystem>
#include <iostream>

#include "cli.hpp"
#include "grid.hpp"
#include "morse_potential.hpp"
#include "morse_static.hpp"
#include "paths.hpp"
#include "taylor.hpp"

namespace fs = std::filesystem;

int main(int argc, char** argv) {
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
    if (!P.quiet) {
        std::cout << "# log saved to: " << csv_path.string() << "\n";
    }

    return 0;
}

