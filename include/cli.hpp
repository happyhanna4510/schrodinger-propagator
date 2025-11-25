#pragma once
#include <string>

struct Params {
    int    N        = 2001;
    double xmax     = 20.0;
    double gamma    = 10.0;
    double Umax     = 0.1;
    bool   Umax_specified = false; //track whether scaling cap was explicitly requested

    std::string init      = "complex-gauss";
    double x0              = 0.0;
    double sigma           = 1.0;
    double k0              = 10.0;
    
    double U0              = 0.0;   // Uniform field amplitude U0=0 (H = H0 + (U0/xmax)*X)
    bool   k0_specified    = false;

    int    eigen_index          = 0;  // [ TEST] default Morse eigenstate index for --init eigen
    bool   eigen_index_specified = false; // [ TEST] track whether user provided --eigen_index




    // evolution parameters
    bool   do_evolve     = false;
    std::string evolve_method = "taylor";
    int    K         = 4;
    double dt        = 1e-5;
    double tmax      = 10.0;
    double tol       = 1e-12;
    int    log_every  = 1000;
    int    csv_every  = 1;
    bool   aggregate  = false;
    int    flush_every = 1000;
    bool   no_theta    = false;
    bool   log_energy  = false;

    int    first     = 10;
    std::string outdir   = "results";
    std::string csv_name;  


    bool evolve_only = false;   // bez generowania potencjału i stanu początkowego
    bool quiet       = false;   // nie wypisywać informacji na konsolę
    bool wide     = false;   // zapisywać |psi|^2 do pliku wide CSV
    bool wide_re     = false;   //  pisać Re[psi] w wide CSV
    bool wide_im     = false;   //  pisać Im[psi] w wide CSV
    bool profile     = false;   // wypisywać profilowanie integratora

    bool export_ref_density = false; // zapisywać gęstość numeryczną i referencyjną do CSV

    std::string cheb_beta_log; // opcjonalny zapis współczynników Czebyszewa

};


Params parse_args(int argc, char** argv);
