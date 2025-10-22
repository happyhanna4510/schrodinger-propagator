#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <iomanip>
#include <vector>
#include <filesystem>

#include "grid.hpp"
#include "morse_potential.hpp"
#include "hamiltonian.hpp"
#include "morse_analytic.hpp"
#include "io.hpp"
#include "taylor.hpp"
#include "snapshots.hpp"

#include "cli.hpp"
#include "paths.hpp"
#include "initial.hpp"
#include "morse_static.hpp"

namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    // 1) Параметры
    Params P = parse_args(argc, argv);

    // 2) Сетка и ИСТИННЫЙ Морс
    Grid g(P.N, P.xmax);
    auto U_true = morse_potential(g, P.gamma);

    // 3) Выходная папка (создаём всегда)
    fs::path out = resolve_outdir(fs::path(argv[0]), P.outdir);
    fs::create_directories(out);

    // 4) Статика Морса: ВСЕГДА пишем (если не попросили только эволюцию)
    if (!P.evolve_only) {
        write_morse_outputs(g, U_true, P.gamma, P.first, out);
    } else {
        std::cout << "# --evolve_only: skipping Morse static outputs.\n";
    }

/////////////////////////////////////
    if (!P.do_evolve) {
    std::cout << "# No evolution (use --evolve [taylor])\n";
    } else {
        auto csv_path = run_taylor_evolution(g, U_true, P, out);
        std::cout << "# log saved to: " << csv_path.string() << "\n";
    }






    /*if (!P.do_evolve) {
    std::cout << "# No evolution (use --evolve [taylor])\n";
    } else {
    // ==== ЭВОЛЮЦИЯ ПО ТЕЙЛОРУ (из конспекта): Umax=0.1, K=4, dt~1e-6 ====

    // 2.1) Масштабируем потенциал к целевому Umax (по конспекту)
    std::vector<double> U_evol = U_true;
    double cur_max = *std::max_element(U_evol.begin(), U_evol.end());
    double scale = (cur_max > 0.0) ? (P.Umax / cur_max) : 1.0;  // P.Umax задай в CLI: --Umax 0.1
    for (auto& v : U_evol) v *= scale;

    std::cout << "# Evolution scaling: current max(U)=" << cur_max
              << ", target Umax=" << P.Umax
              << ", scale=" << scale << "\n";

    // 2.2) Гамильтониан для эволюции и трёхдиагональ
    Eigen::MatrixXd H_evol = build_hamiltonian(g, U_evol);

    // если у тебя есть оптимальный трёхдиагональный формат — используем его:
    Tridiag T = make_tridiag_from_dense(H_evol);  // твоя функция

    // 2.3) Начальная волна: ГАУСС (выбери один вариант)

    // Вариант А: стоящий гаусс у минимума (x0~1.0, sigma~0.35)
    Eigen::VectorXcd psi_init = gaussian_on_inner(g, 1.0, 0.35).cast<std::complex<double>>();

    // Вариант Б: движущийся гаусс (раскомментируй, если хочешь движение)
    // Eigen::VectorXcd psi_init = gaussian_complex_on_inner(g, 1.0, 0.35, 1.0); // k0=1.0

    // 2.4) Параметры Тейлора (конспект: K=4, маленький dt)
    int    K      = (P.K > 0 ? P.K : 4);           // или CLI: --K 4
    double dt     = (P.dt > 0 ? P.dt : 1e-6);      // --dt 1e-6
    double tmax   = (P.tmax > 0 ? P.tmax : 1e-3);  // --tmax 1e-3 (для старта)
    int    nsteps = (int)std::llround(tmax / dt);

    // лог и снапшоты
    int log_every = (P.log_every > 0 ? P.log_every : 100);
    fs::path csv_path = out / (P.csv_name.empty() ? "gauss_taylor.csv" : P.csv_name);

    // ось X для вывода
    std::vector<double> x_inner(g.N-2); for (int i=0;i<g.N-2;++i) x_inner[i] = g.x[i+1];
    const std::vector<double>* x_ptr = &x_inner;

    // 2.5) evolve_taylor_tridiag — твоя функция; в старом коде она принимала phi0 и E0.
    Eigen::VectorXcd phi0 = psi_init;  // сохраним исходное для контроля
    double E0 = 0.0;                   // для общего пакета не нужен, но если сигнатура требует — передаём 0.

    std::cout << "\n# === Time evolution (Taylor K="<<K
              << ", dt="<<dt<<", tmax="<<tmax<<") ===\n";

    evolve_taylor_tridiag(
        T,                // трёхдиагональный гамильтониан
        psi_init,         // текущее psi (внутри будет обновляться)
        phi0,             // копия начального состояния (для метрик/логов)
        E0,               // не используется для общего гаусса
        g.dx,             // шаг сетки для нормы
        dt, nsteps, K,    // время и порядок ряда
        log_every,        // раз в сколько шагов лог/проверка нормы
        csv_path.string(),// файл лога
        x_ptr,            // передаём X, чтобы снепшоты писались (если поддерживается)
        ""                // snaps_dir — пусто, если не нужен покадровый вывод
    );

    std::cout << "# log saved to: " << csv_path.string() << "\n";
}
*/

    return 0;
}
