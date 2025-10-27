#include "cli.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace {

void warn_missing(const char* flag) {
    std::cerr << "warning: missing value after " << flag << "\n";
}

}

Params parse_args(int argc, char** argv) {
    Params p;

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

        if      (s == "--N")        geti(i, p.N);
        else if (s == "--xmax")     getd(i, p.xmax);
        else if (s == "--gamma")    getd(i, p.gamma);

        else if (s == "--Umax")     getd(i, p.Umax);
        else if (s == "--Vcap")     getd(i, p.Umax); // alias

        else if (s == "--first")    geti(i, p.first);

        else if (s == "--evolve") { p.do_evolve = true; gets(i, p.evolve_method); }
        else if (s == "--dt")       getd(i, p.dt);
        else if (s == "--tmax")     getd(i, p.tmax);
        else if (s == "--K")        geti(i, p.K);
        else if (s == "--tol")      getd(i, p.tol);
        else if (s == "--log" || s == "--log-every") geti(i, p.log_every);
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

        else {
            std::cerr << "warning: unknown option: " << s << "\n";
        }
    }

    return p;
}

