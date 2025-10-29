#pragma once
#include <string>

struct Params {
    int    N        = 2001;
    double xmax     = 30.0;
    double gamma    = 10.0;
    double Umax     = 0.1;

    // evolution parameters
    bool   do_evolve     = false;
    std::string evolve_method = "taylor";
    int    K         = 4;
    double dt        = 1e-5;
    double tmax      = 10.0;
    double tol       = 1e-12;
    int    log_every  = 10000;
    int    csv_every  = 1;
    bool   aggregate  = false;
    int    flush_every = 1000;
    bool   no_theta    = false;

    int    first     = 10;
    std::string outdir   = "results";
    std::string csv_name;  


    bool evolve_only = false;   // bez generowania potencjału i stanu początkowego
    bool quiet       = false;   // nie wypisywać informacji na konsolę
    bool wide     = false;   // zapisywać |psi|^2 do pliku wide CSV
    bool wide_re     = false;   //  pisać Re[psi] w wide CSV
    bool wide_im     = false;   //  pisać Im[psi] w wide CSV

};


Params parse_args(int argc, char** argv);
