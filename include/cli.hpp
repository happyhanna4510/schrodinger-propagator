#pragma once
#include <string>

struct Params {
    int    N        = 2001;
    double xmax     = 30.0;
    double gamma    = 10.0;
    double Umax     = 0.1;

    // эволюция
    bool   do_evolve     = false;
    std::string evolve_method = "taylor";
    int    K         = 4;
    double dt        = 1e-5;
    double tmax      = 10.0;
    int    log_every = 10000;

    // выводы
    int    first     = 10;
    std::string outdir   = "results";
    std::string csv_name;  // пусто по умолчанию


    // НОВЫЕ ФЛАГИ
    bool evolve_only = false;   // пропустить статическую часть, если морс уже есть
    bool quiet       = false;   // не печатать таблицы/мониторинг
    bool no_wide     = false;   // не писать psi_*_wide.csv

};


Params parse_args(int argc, char** argv);
