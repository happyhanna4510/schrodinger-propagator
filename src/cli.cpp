#include "cli.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace {

void warn_missing(const char* flag) {
    std::cerr << "warning: missing value after " << flag << "\n";
}

void print_help(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n\n"
              << "Initial state options:\n"
              << "  --init <complex-gauss|real-gauss>  Select Gaussian type (default complex-gauss).\n"
              << "  --x0 <value>                        Gaussian center (default 0).\n"
              << "  --sigma <value>                     Gaussian width (default 1). Must be > 0.\n"
              << "  --k0 <value>                        Plane-wave wavenumber for complex-gauss (default 10).\n"
              << "  --U0 <value>                        Uniform field amplitude in H = H0 + (U0/xmax) * X (default 0).\n\n"
              << "Evolution options:\n"
              << "  --evolve <method>                   Enable time evolution using taylor|rk4|cheb.\n"
              << "  --dt <value>                        Time step size.\n"
              << "  --tmax <value>                      Total simulated time.\n"
              << "  --K <value>                         Taylor/Chebyshev order limit.\n"
              << "  --tol <value>                       Chebyshev tolerance.\n\n"
              << "Grid & potential:\n"
              << "  --N <value>                         Number of grid points.\n"
              << "  --xmax <value>                      Half-width of domain [-xmax, xmax].\n"
              << "  --gamma <value>                     Morse potential parameter.\n"
              << "  --Umax <value>                      Target max potential amplitude for scaling.\n\n"
              << "Output control:\n"
              << "  --outdir <dir>                      Output directory (default results).\n"
              << "  --csv <name>                        Custom CSV filename.\n"
              << "  --log-every <steps>                 Console log cadence.\n"
              << "  --log-energy                       Enable energy logging (writes separate CSV).\n"
              << "  --csv-every <steps>                 CSV log cadence.\n"
              << "  --flush-every <rows>                CSV flush cadence.\n"
              << "  --aggregate                         Aggregate log timing stats.\n"
              << "  --no-theta                          Skip theta metrics.\n"
              << "  --wide[,-re,-im]                    Emit wide output dumps.\n"
              << "  --export-ref-density                Save numerical and reference densities to CSV.\n"
              << "  --quiet                             Suppress console logs.\n"
              << "  --evolve_only                       Skip regenerating static Morse outputs.\n"
              << "  --profile                           Print integrator profiling data.\n"
              << "  --help                              Show this message and exit.\n\n"
              << "Examples:\n"
              << "  " << prog << " --evolve taylor --init complex-gauss --x0 -10 --sigma 1.5 --k0 12\n"
              << "  " << prog << " --evolve taylor --init complex-gauss --x0 -10 --sigma 1.2 --k0 15 --U0 2.0\n"
              << "  " << prog << " --evolve rk4 --init real-gauss --x0 0 --sigma 2.0\n";
}

std::string canonicalize_init(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    std::replace(value.begin(), value.end(), '_', '-');
    return value;
}

}

Params parse_args(int argc, char** argv) {
    Params p;
    bool k0_specified = false;

    auto getd = [&](int& i, double& dst) {
        if (i + 1 < argc) dst = std::atof(argv[++i]);
        else warn_missing(argv[i]);
    };
    auto geti = [&](int& i, int& dst) {
        if (i + 1 < argc) dst = std::atoi(argv[++i]);
        else warn_missing(argv[i]);
    };
    auto gets = [&](int& i, std::string& s) {
        if (i + 1 < argc) s = argv[++i];
        else warn_missing(argv[i]);
    };

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];

        if (s == "--help" || s == "-h") {
            print_help(argv[0]);
            std::exit(0);
        }

        if      (s == "--N")        geti(i, p.N);
        else if (s == "--xmax")     getd(i, p.xmax);
        else if (s == "--gamma")    getd(i, p.gamma);

        else if (s == "--Umax")  { getd(i, p.Umax); p.Umax_specified = true; }
        else if (s == "--Vcap")  { getd(i, p.Umax); p.Umax_specified = true; } // alias

        else if (s == "--init") { gets(i, p.init); }
        else if (s == "--eigen_index" || s == "--eigen-index") { geti(i, p.eigen_index); p.eigen_index_specified = true; } // [ TEST] parse eigen index selector

        else if (s == "--x0")       getd(i, p.x0);
        else if (s == "--sigma")    getd(i, p.sigma);
        else if (s == "--k0")      { getd(i, p.k0); k0_specified = true; }
        else if (s == "--U0")       getd(i, p.U0);

        else if (s == "--first")    geti(i, p.first);

        else if (s == "--evolve") { p.do_evolve = true; gets(i, p.evolve_method); }
        else if (s == "--dt")       getd(i, p.dt);
        else if (s == "--tmax")     getd(i, p.tmax);
        else if (s == "--K")        geti(i, p.K);
        else if (s == "--tol")      getd(i, p.tol);
        else if (s == "--log" || s == "--log-every") geti(i, p.log_every);
        else if (s == "--log-energy") p.log_energy = true;
        else if (s == "--csv-every") geti(i, p.csv_every);
        else if (s == "--flush-every") geti(i, p.flush_every);

        else if (s == "--csv") {
            if (i + 1 < argc) p.csv_name = argv[++i];
            else warn_missing(argv[i]);
        }
        else if (s == "--outdir")   gets(i, p.outdir);

        else if (s == "--evolve_only" || s == "--evolve-only") p.evolve_only = true;
        else if (s == "--quiet")       p.quiet       = true;
        else if (s == "--wide" ) p.wide = true;
        //else if (s == "--wide-re")     p.wide_re     = true;
        //else if (s == "--wide-im")     p.wide_im     = true;

        else if (s == "--aggregate") p.aggregate = true;
        else if (s == "--no-theta") p.no_theta = true;
        else if (s == "--profile") p.profile = true;
        else if (s == "--export-ref-density") p.export_ref_density = true;

        else {
            std::cerr << "warning: unknown option: " << s << "\n";
        }
    }

    p.init = canonicalize_init(p.init);
    if (p.init != "complex-gauss" && p.init != "real-gauss" && p.init != "eigen") {
        std::cerr << "error: --init must be 'complex-gauss', 'real-gauss', or 'eigen' (got '" << p.init << "')\n"; // [AI PATCH TEST] allow eigen init option
        std::exit(1);
    }

    if (p.init == "eigen" && !p.eigen_index_specified) {
        p.eigen_index = 0; // [ TEST] default to ground state when user omits --eigen_index
    }

    if (p.sigma <= 0.0) {
        std::cerr << "error: sigma must be > 0 (got " << p.sigma << ")\n";
        std::exit(1);
    }

    p.k0_specified = k0_specified;

    return p;

}

