# plot_morse.py
import argparse, csv, os, math
import numpy as np
import matplotlib.pyplot as plt

def read_xy_csv(path):
    xs, ys = [], []
    with open(path, "r", newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        for row in r:
            if not row: 
                continue
            try:
                xs.append(float(row[0]))
                ys.append(float(row[1]))
            except ValueError:
                continue
    return np.array(xs), np.array(ys)

def morse_nmax(gamma: float) -> int:
    return max(int(math.floor(gamma - 0.5 - 1e-12)), -1)

def morse_energies(gamma: float, how_many: int):
    nmax = morse_nmax(gamma)
    m = min(how_many, max(nmax + 1, 0))
    E = []
    for n in range(m):
        En = 1.0 - ((gamma - (n + 0.5))**2) / (gamma*gamma)
        E.append(En)
    return np.array(E)

def autodomain(x, y, xleft=None, xright=None, xmax_sym=None, ymin=None, ymax=None, extra_y=None):
    # X-limits
    if xleft is not None and xright is not None:
        xl, xr = xleft, xright
    elif xmax_sym is not None:
        xl, xr = -abs(xmax_sym), abs(xmax_sym)
    else:
        pad = 0.05*(x.max()-x.min())
        xl, xr = x.min()-pad, x.max()+pad

    # ограничим массив по выбранному X, чтобы подобрать Y по видимой части
    mask = (x >= xl) & (x <= xr)
    yvis = y[mask] if mask.any() else y

    # Y-limits
    if ymin is None or ymax is None:
        cand = [yvis.min(), yvis.max()]
        if extra_y is not None and len(extra_y):
            cand += [extra_y.min(), extra_y.max()]
        y_min, y_max = min(cand), max(cand)
        if y_min == y_max:
            y_min -= 1.0; y_max += 1.0
        pad = 0.05*(y_max - y_min)
        ylo = y_min - pad if ymin is None else ymin
        yhi = y_max + pad if ymax is None else ymax
    else:
        ylo, yhi = ymin, ymax

    return (xl, xr), (ylo, yhi)

def main():
    p = argparse.ArgumentParser(description="Plot Morse potential with optional energy levels.")
    p.add_argument("--csv", default=None, help="путь к morse_potential.csv")
    # диапозоны осей
    p.add_argument("--xmax", type=float, default=None,
                   help="симметричный диапазон по X: покажет от -xmax до +xmax")
    p.add_argument("--xleft", type=float, default=None, help="левая граница X (если задана, перекроет --xmax)")
    p.add_argument("--xright", type=float, default=None, help="правая граница X (если задана, перекроет --xmax)")
    p.add_argument("--ymin", type=float, default=None)
    p.add_argument("--ymax", type=float, default=None)
    # уровни
    p.add_argument("--gamma", type=float, default=None, help="если задано — нарисовать уровни энергии")
    p.add_argument("--levels", type=int, default=None, help="сколько уровней рисовать (по умолчанию максимум)")
    p.add_argument("--no-levels", action="store_true", help="не рисовать уровни даже если задана gamma")

    
    p.add_argument("--asymptote", action="store_true", help="показать горизонтальную асимптоту")

    # оформление
    p.add_argument("--figw", type=float, default=7.0)
    p.add_argument("--figh", type=float, default=4.2)
    p.add_argument("--dpi", type=int, default=140)
    p.add_argument("--grid", action="store_true")
    p.add_argument("--title", default="Morse potential")
    p.add_argument("--save", default=None, help="PNG путь, напр. figures/morse_g6.png")
    args = p.parse_args()

    # автопоиск CSV
    candidates = []
    if args.csv: candidates.append(args.csv)
    candidates += [
        os.path.join("results", "morse_potential.csv"),
        os.path.join("..", "results", "morse_potential.csv"),
        os.path.join("build","Release","morse_potential.csv"),
        os.path.join("build","morse_potential.csv"),
    ]
    csv_path = next((c for c in candidates if os.path.exists(c)), None)
    if csv_path is None:
        raise SystemExit("Не найден morse_potential.csv (сначала запусти exe с --outdir results)")

    x, U = read_xy_csv(csv_path)

    # уровни энергии
    Es = None
    if not args.no_levels and args.gamma is not None:
        nmax = morse_nmax(args.gamma)
        L = (nmax + 1) if args.levels is None else args.levels
        Es = morse_energies(args.gamma, L)

    # авто-пределы осей
    (xl, xr), (yl, yr) = autodomain(
        x, U,
        xleft=args.xleft, xright=args.xright, xmax_sym=args.xmax,
        ymin=args.ymin, ymax=args.ymax,
        extra_y=Es
    )

    # рисуем
    plt.figure(figsize=(args.figw, args.figh), dpi=args.dpi)
    plt.plot(x, U, label="Morse U(x)")
    if args.grid:
        plt.grid(True, alpha=0.25)
    plt.axhline(0, lw=1, ls="--")

    # горизонтальная асимптота
    if args.asymptote:
        y_as = 1.0          # безразмерная форма: U(+∞)=1
        # если у тебя есть физический режим с De:
        # if args.phys: y_as = args.De
        plt.axhline(y_as, lw=1.2, ls="--", alpha=0.8)
        plt.text(xl + 0.02*(xr-xl), y_as + 0.02*(yr-yl),
                f"asymptote = {y_as:g}", va="bottom")


    # уровни
    if Es is not None and len(Es):
        for E in Es:
            plt.hlines(E, xl, xr, linewidth=1)
        # подпись над верхним уровнем
        xspan = xr - xl
        yspan = yr - yl
        plt.text(xl + 0.05*xspan, Es[-1] + 0.03*yspan,
                 f"levels: {len(Es)} (gamma={args.gamma})")

    plt.xlim(xl, xr)
    plt.ylim(yl, yr)
    plt.title(args.title)
    plt.xlabel("x")
    plt.ylabel("U(x)")
    plt.legend(loc="best")
    plt.tight_layout()

    if args.save:
        os.makedirs(os.path.dirname(args.save), exist_ok=True)
        plt.savefig(args.save, dpi=max(200, args.dpi))
    plt.show()

if __name__ == "__main__":
    main()
