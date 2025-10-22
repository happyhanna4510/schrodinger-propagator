#include "cli.hpp"

#include <cstdlib>    // std::atof, std::atoi
#include <cstring>    // strcmp (если используешь)
#include <string>     // std::string
#include <iostream>   // std::cerr, std::cout

Params parse_args(int argc, char** argv){
    Params p;

    auto getd = [&](int& i, double& dst){
        if (i+1 < argc) dst = std::atof(argv[++i]);
        else std::cerr << "warning: missing value after " << argv[i] << "\n";
    };
    auto geti = [&](int& i, int& dst){
        if (i+1 < argc) dst = std::atoi(argv[++i]);
        else std::cerr << "warning: missing value after " << argv[i] << "\n";
    };
    auto gets = [&](int& i, std::string& s){
        if (i+1 < argc) s = argv[++i];
        else std::cerr << "warning: missing value after " << argv[i] << "\n";
    };

    for (int i = 1; i < argc; ++i){
        std::string s = argv[i];

        // базовые параметры сетки/потенциала
        if      (s == "--N")        geti(i, p.N);
        else if (s == "--xmax")     getd(i, p.xmax);
        else if (s == "--gamma")    getd(i, p.gamma);

        // cap/scale параметр (сохраняю оба алиаса, если у вас в структуре есть p.Umax)
        else if (s == "--Umax")     getd(i, p.Umax);
        else if (s == "--Vcap")     getd(i, p.Umax);     // [alias] чтобы старые скрипты не сломались

        // сколько уровней печатать/сохранять
        else if (s == "--first")    geti(i, p.first);

        // эволюция
        else if (s == "--evolve") { p.do_evolve = true; gets(i, p.evolve_method); }
        else if (s == "--dt")       getd(i, p.dt);
        else if (s == "--tmax")     getd(i, p.tmax);
        else if (s == "--K")        geti(i, p.K);
        else if (s == "--log")      geti(i, p.log_every);

        // вывод путей
        else if (s == "--csv" && i+1 < argc) p.csv_name = argv[++i];
        else if (s == "--outdir")   gets(i, p.outdir);

        // флаги режимов
        else if (s == "--evolve_only" || s == "--evolve-only") p.evolve_only = true;  // пропустить статику если morse_* уже есть
        else if (s == "--quiet")       p.quiet       = true;  // не печатать таблицы/прогресс
        else if (s == "--no_wide" || s == "--no-wide") p.no_wide = true;  // не писать psi_*_wide.csv

        // неизвестные ключи — просто предупредим
        else {
            std::cerr << "warning: unknown option: " << s << "\n";
        }
    }

    return p;
}
