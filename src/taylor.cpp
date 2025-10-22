#include "taylor.hpp"
#include <fstream>
#include <cstdio>
#include "snapshots.hpp"
#include <filesystem>
#include <iomanip> // для setprecision

#include "cli.hpp"
#include "hamiltonian.hpp"
#include "initial.hpp"      // gaussian_on_inner / gaussian_complex_on_inner
#include "io.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <Eigen/Core>
#include <sstream>  // вверху файла, если ещё не подключён


namespace fs = std::filesystem;

static const cplx I(0.0, 1.0);

inline double l2_norm_dx(const Eigen::VectorXcd& v, double dx){
    long double s=0.0L; for(Eigen::Index i=0;i<v.size();++i) s+=std::norm(v[i]);
    return (double)(s * (long double)dx);
}
inline std::complex<double> inner_dx(const Eigen::VectorXcd& a, const Eigen::VectorXcd& b, double dx){
    return dx * a.dot(b);
}

TaylorLogRow metrics_vs_ground(const Eigen::VectorXcd& psi,
                               const Eigen::VectorXcd& phi0,
                               double E0, double t, double dx, int step)
{
    const double npsi  = l2_norm_dx(psi, dx);
    const double nphi0 = l2_norm_dx(phi0, dx);
    const auto   c     = inner_dx(phi0, psi, dx);
    const double p0    = (npsi>0 && nphi0>0) ? (std::norm(c)/(npsi*nphi0)) : 0.0;

    const auto psi_exact = std::exp(std::complex<double>(0,-1)*E0*t) * phi0; // эталон для φ0
    const double err = (psi - psi_exact).norm(); // ок оставить евклидову

    return {step, t, npsi, p0, err};
}


// ===== трёхдиагональная обвязка =====
Tridiag make_tridiag_from_dense(const Eigen::MatrixXd& H){
    const int M = (int)H.rows();
    Tridiag T{Eigen::VectorXd(M), Eigen::VectorXd(M-1), Eigen::VectorXd(M-1)};
    T.a = H.diagonal();
    T.b = H.diagonal(1);
    T.c = H.diagonal(-1);
    return T;
}

// ===== восстановление плотной матрицы из трёхдиагональной =====
Eigen::MatrixXd tridiag_to_dense(const Tridiag& T) {
    const int M = static_cast<int>(T.a.size());
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(M, M);
    H.diagonal()      = T.a;
    H.diagonal(1)     = T.b;
    H.diagonal(-1)    = T.c;
    return H;
}

// y = H*x для трёхдиагональной H (вещественные диагонали; x,y — комплексные)
static inline void tridiag_mul(const Tridiag& T,
                               const Eigen::VectorXcd& x,
                               Eigen::VectorXcd& y)
{
    const int M = (int)T.a.size();
    y.setZero(M);

    if (M == 0) return;
    if (M == 1) { y[0] = T.a[0]*x[0]; return; }

    y[0] = T.a[0]*x[0] + T.b[0]*x[1];
    for (int i=1; i<M-1; ++i){
        y[i] = T.c[i-1]*x[i-1] + T.a[i]*x[i] + T.b[i]*x[i+1];
    }
    y[M-1] = T.c[M-2]*x[M-2] + T.a[M-1]*x[M-1];
}

// один шаг Тейлора порядка K через трёхдиагональное умножение
void taylor_step_tridiag(const Tridiag& T, Eigen::VectorXcd& psi, double dt, int K){
    Eigen::VectorXcd sum = psi;
    Eigen::VectorXcd vk  = psi;
    Eigen::VectorXcd tmp(psi.size());

    for (int k=1; k<=K; ++k){
        tridiag_mul(T, vk, tmp);              // tmp = H * vk
        vk = (-I * dt / double(k)) * tmp;     // vk = (-i dt / k) * tmp
        sum += vk;
    }
    psi.swap(sum);
}


// находит основное состояние для H(g, U) и нормирует собственные векторы по dx
static std::pair<double, Eigen::VectorXcd>
ground_state_for(const Grid& g, const std::vector<double>& U)
{
    Eigen::MatrixXd H = build_hamiltonian(g, U);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    if (es.info()!=Eigen::Success) {
        throw std::runtime_error("eigensolve failed");
    }
    Eigen::VectorXd Evals = es.eigenvalues();
    Eigen::MatrixXd EV    = es.eigenvectors();
    renormalize_eigenvectors(EV, g.dx);                  // ВАЖНО: учёт dx!

    int idx0 = 0;                                        // минимум энергии
    Eigen::VectorXcd phi0 = EV.col(idx0).cast<std::complex<double>>();
    double E0 = Evals(idx0);

    // фиктивная фаза, чтобы не прыгал знак
    if (phi0(0).real() < 0) phi0 = -phi0;

    return {E0, phi0};
}


static inline double l2norm(const Eigen::VectorXcd& v, double dx){
    // ||v||_L2 ~ sqrt( sum |v_i|^2 * dx )
    return std::sqrt( v.array().abs2().sum() * dx );  // скобки после array()
}

RevTestResult test_time_reversibility(const Tridiag& T,
                                      const Eigen::VectorXcd& psi_init,
                                      double dx, double dt, int nsteps, int K)
{
    RevTestResult r{};
    r.norm0 = l2norm(psi_init, dx);

    Eigen::VectorXcd psi = psi_init;

    auto upd_drift = [&](double& cur_max, const Eigen::VectorXcd& v){
        double n = l2norm(v, dx);
        cur_max = std::max(cur_max, std::abs(n - 1.0));
    };

    // --- прямой прогон +dt ---
    double max_drift = 0.0;
    for (int i = 0; i < nsteps; ++i){
        taylor_step_tridiag(T, psi, dt, K);
        if ((i & 0x3FF) == 0) upd_drift(max_drift, psi); // раз в 1024 шага
    }
    r.norm_fwd = l2norm(psi, dx);
    upd_drift(max_drift, psi);

    // --- обратный прогон -dt ---
    for (int i = 0; i < nsteps; ++i){
        taylor_step_tridiag(T, psi, -dt, K);
        if ((i & 0x3FF) == 0) upd_drift(max_drift, psi);
    }
    r.norm_back = l2norm(psi, dx);
    upd_drift(max_drift, psi);

    // ошибка возврата
    r.err_back = (psi - psi_init).norm(); // евклидова (без dx) — ок для сравнения
    r.max_drift = max_drift;
    return r;
}


static double prob_slice(const Eigen::VectorXcd& psi, int i0, int i1, double dx){
    double s = 0.0;
    for (int i=i0; i<i1; ++i) s += std::norm(psi(i));
    return s * dx;
}


// эволюция + лог в CSV
void evolve_taylor_tridiag(const Tridiag& T,
                           const Eigen::VectorXcd& psi_init,
                           const Eigen::VectorXcd& phi0,
                           double E0, double dx, double dt, int nsteps, int K,
                           int log_every, const std::string& csv_path,
                           const std::vector<double>* x_inner,
                           const std::string& /*snaps_dir - не нужен теперь*/,
                           LogExtras extras)    // <=== НОВОЕ
{
    std::ofstream f;
    if (!csv_path.empty()){
        f.open(csv_path, std::ios::out);

        // Заголовок CSV — без p0 и err_phi0 по умолчанию
        // при необходимости можно вернуть их, включив флаги
        std::ostringstream h;
        h << "step,t,norm";
        if (has(extras, LogExtras::P0))     h << ",p0";
        if (has(extras, LogExtras::ErrPhi)) h << ",err_phi0";
        h << ",err_true,theta,Pleft,Pright\n";  // <— добавили
        f << h.str();
    }

    io::WideDump wide(csv_path, x_inner /*, false, false*/); // по умолчанию только abs2


    Eigen::VectorXcd psi = psi_init;
    double t = 0.0;

    // === подготовка точного базиса и начальных коэф. c0 ===
    Eigen::MatrixXd H = tridiag_to_dense(T);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    Eigen::VectorXd E = es.eigenvalues();
    Eigen::MatrixXd V = es.eigenvectors();
    renormalize_eigenvectors(V, dx);

    Eigen::VectorXcd c0 = (dx * V.transpose()) * psi_init;

    double Emin = E.minCoeff(), Emax = E.maxCoeff(), dE = Emax - Emin;
    std::cout << "# Spectrum: Emin=" << Emin << "  Emax=" << Emax
              << "  dE=" << dE << "\n";

    static double norm0 = -1.0;

    for (int step = 0; step <= nsteps; ++step)
    {
        if (step % log_every == 0)
        {
            auto m = metrics_vs_ground(psi, phi0, E0, t, dx, step); // m.p0 и m.err_exact доступны

            if (norm0 < 0.0) norm0 = m.norm;
            double drift = m.norm - norm0;

            // общая фаза
            Eigen::ArrayXcd phase = (-std::complex<double>(0,1) * E.array() * t).exp();

            // точная волна из c0
            Eigen::VectorXcd psi_exact =
                V.cast<std::complex<double>>() * (c0.array() * phase).matrix();
            double err_true = (psi - psi_exact).norm();

            // коэффициенты из текущей волны и точные коэф.
            Eigen::VectorXcd c_num   = (dx * V.transpose()) * psi;
            Eigen::VectorXcd c_exact = (c0.array() * phase).matrix();

            // Θ(t) = sum |c_exact - c_num|^2
            double theta = (c_exact - c_num).squaredNorm();

            // ---------- КОНСОЛЬ: без p0/err_phi0 по умолчанию ----------
            std::ostringstream line;
            line << std::setw(8) << m.step
                 << "  t=" << std::scientific << std::setprecision(6) << m.t
                 << "  norm=" << std::setprecision(15) << m.norm
                 << "  drift=" << std::scientific << std::setprecision(3) << drift;

            if (has(extras, LogExtras::P0))
                line << "  p0=" << std::fixed << std::setprecision(10) << m.p0;

            if (has(extras, LogExtras::ErrPhi))
                line << "  err_phi0=" << std::scientific << std::setprecision(3) << m.err_exact;

            line << "  err_true=" << std::scientific << std::setprecision(3) << err_true
                 << "  Theta="    << std::scientific << std::setprecision(3) << theta;

            std::cout << line.str() << "\n";


            int M = psi.size();
            int k = 10; // ширина буфера (можешь менять)
            double P_left  = prob_slice(psi, 0, k, dx);
            double P_right = prob_slice(psi, M-k, M, dx);


            // ---------- CSV: те же поля, что и в заголовке ----------
            if (f.is_open())
            {
                std::ostringstream row;
                row << m.step << ','
                    << std::setprecision(15) << m.t << ','
                    << m.norm;

                if (has(extras, LogExtras::P0))
                    row << ',' << std::setprecision(10) << m.p0;

                if (has(extras, LogExtras::ErrPhi))
                    row << ',' << std::scientific << std::setprecision(3) << m.err_exact;

                row << ',' << std::scientific << std::setprecision(3) << err_true
                    << ',' << std::scientific << std::setprecision(3) << theta
                    << ',' << std::scientific << std::setprecision(6) << P_left
                    << ',' << std::scientific << std::setprecision(6) << P_right
                    << '\n';

                f << row.str();
            }

            wide.write(psi, t);
            
            std::cout << "  edge(P): L=" << std::scientific << std::setprecision(3) << P_left
          << " R=" << P_right << "\n";

        }

        if (step == nsteps) break;
        taylor_step_tridiag(T, psi, dt, K);
        t += dt;
    }
}



///////////////////////////////////////////////////////////////////////////
fs::path run_taylor_evolution(const Grid& g,
                              const std::vector<double>& U_true,
                              const Params& P,
                              const fs::path& out_dir)
{
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
    int    log_every = (P.log_every > 0 ? P.log_every : 100);


    // 5) Имя CSV: если задано --csv, используем его, иначе собираем автоматически
    fs::create_directories(out_dir);
    // 5) Имя CSV: если задано --csv, используем его, иначе авто
    fs::path csv_path = io::make_csv_path(out_dir, P, K, dt, g);


    // ось X для вывода
    std::vector<double> x_inner(g.N-2); for (int i=0;i<g.N-2;++i) x_inner[i] = g.x[i+1];
    const std::vector<double>* x_ptr = &x_inner;

    // 2.5) evolve_taylor_tridiag — твоя функция; в старом коде она принимала phi0 и E0.
    auto [E0, phi0] = ground_state_for(g, U_evol);


    // === тест обратимости ===
{
    RevTestResult rt = test_time_reversibility(T, psi_init, g.dx, dt, nsteps, K);
    std::cout << "\n# Reversibility test (K=" << K << ", dt=" << dt
              << ", nsteps=" << nsteps << ")\n"
              << "  err_back = " << std::scientific << rt.err_back << "\n"
              << "  norm0    = " << std::setprecision(12) << rt.norm0 << "\n"
              << "  norm_fwd = " << std::setprecision(12) << rt.norm_fwd << "\n"
              << "  norm_back= " << std::setprecision(12) << rt.norm_back << "\n"
              << "  max|norm-1| (both passes) = " << std::scientific << rt.max_drift << "\n";
}
// === конец теста ===


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
        ""  ,              // snaps_dir — пусто, если не нужен покадровый вывод
        LogExtras::None
    );
    
    return csv_path;
}

