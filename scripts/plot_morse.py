# plot_morse.py
import argparse, csv, os, math
import numpy as np
import matplotlib.pyplot as plt

def read_xy_csv(path):
    # Czytanie prostego CSV z dwiema kolumnami: x, y (pomijamy puste/nieparsowalne wiersze)
    xs, ys = [], []
    with open(path, "r", newline="", encoding="utf-8") as f:
        r = csv.reader(f)
        for row in r:
            if not row:
                continue
            try:
                xs.append(float(row[0])); ys.append(float(row[1]))
            except ValueError:
                continue
    return np.array(xs), np.array(ys)

def morse_nmax(gamma: float) -> int:
    # Liczba stanów związanych: floor(gamma - 1/2); -1 oznacza „brak”
    return max(int(math.floor(gamma - 0.5 - 1e-12)), -1)

def morse_energies(gamma: float, how_many: int):
    # Energie w postaci bezwymiarowej (U(+∞)=1)
    nmax = morse_nmax(gamma)
    m = min(how_many, max(nmax + 1, 0))
    E = []
    for n in range(m):
        En = 1.0 - ((gamma - (n + 0.5))**2) / (gamma*gamma)
        E.append(En)
    return np.array(E)

def autodomain(x, y, xleft=None, xright=None, xmax_sym=None, ymin=None, ymax=None, extra_y=None):
    # Dobór zakresu X: [xleft, xright] > symetryczny [-xmax_sym, +xmax_sym] > auto (5% margines)
    if xleft is not None and xright is not None:
        xl, xr = xleft, xright
    elif xmax_sym is not None:
        xl, xr = -abs(xmax_sym), abs(xmax_sym)
    else:
        pad = 0.05*(x.max()-x.min()); xl, xr = x.min()-pad, x.max()+pad
    mask = (x >= xl) & (x <= xr)
    yvis = y[mask] if mask.any() else y
    if ymin is None or ymax is None:
        cand = [yvis.min(), yvis.max()]
        if extra_y is not None and len(extra_y):
            cand += [extra_y.min(), extra_y.max()]
        y_min, y_max = min(cand), max(cand)
        if y_min == y_max: y_min -= 1.0; y_max += 1.0
        pad = 0.05*(y_max - y_min)
        ylo = y_min - pad if ymin is None else ymin
        yhi = y_max + pad if ymax is None else ymax
    else:
        ylo, yhi = ymin, ymax
    return (xl, xr), (ylo, yhi)

def main():
    p = argparse.ArgumentParser(description="Wykres potencjału Morse’a z opcjonalnymi poziomami energii.")
    p.add_argument("--csv", default=None, help="Ścieżka do pliku morse_potential.csv")
    # zakresy osi
    p.add_argument("--xmax", type=float, default=None, help="Symetryczny zakres po X: od -xmax do +xmax")
    p.add_argument("--xleft", type=float, default=None, help="Lewa granica X (przeważa nad --xmax)")
    p.add_argument("--xright", type=float, default=None, help="Prawa granica X (przeważa nad --xmax)")
    p.add_argument("--ymin", type=float, default=None, help="Dolna granica Y")
    p.add_argument("--ymax", type=float, default=None, help="Górna granica Y")
    # poziomy
    p.add_argument("--gamma", type=float, default=None, help="Jeśli podane — narysować poziomy energii")
    p.add_argument("--levels", type=int, default=None, help="Ile poziomów (domyślnie maksimum)")
    p.add_argument("--no-levels", action="store_true", help="Nie rysować poziomów nawet jeśli jest --gamma")
    # dodatki
    p.add_argument("--asymptote", action="store_true", help="Pokaż asymptotę U(+∞)=1")
    # wygląd
    p.add_argument("--figw", type=float, default=7.0, help="Szerokość figury [cale]")
    p.add_argument("--figh", type=float, default=4.2, help="Wysokość figury [cale]")
    p.add_argument("--dpi", type=int, default=140, help="DPI figury")
    p.add_argument("--grid", action="store_true", help="Włącz siatkę")
    p.add_argument("--title", default="Potencjał Morse’a", help="Tytuł wykresu")
    p.add_argument("--save", default=None, help="Ścieżka PNG, np. figures/morse_g6.png")
    # root dla layoutu results/morse/g10
    p.add_argument("--root", default=os.path.join("results", "morse"),
                   help="Katalog bazowy z danymi (domyślnie results/morse). Szuka podfolderów g<gamma>.")
    args = p.parse_args()


        # --- Normalize root to handle running from scripts/ ---
    here = os.path.dirname(os.path.abspath(__file__))
    root_candidates = [
        args.root,
        os.path.join(here, args.root),
        os.path.join(here, "..", args.root),
        os.path.join(os.getcwd(), args.root),
    ]
    roots = [os.path.abspath(p) for p in root_candidates if os.path.isdir(os.path.abspath(p))]
    root_norm = roots[0] if roots else os.path.abspath(args.root)


    # --- Autodetekcja CSV (внутри main!) ---
    candidates = []
    if args.csv:
        candidates.append(args.csv)
    if args.gamma is not None:
        g_int = str(int(round(args.gamma)))
        for d in [
            os.path.join(root_norm , f"g{g_int}"),
            os.path.join(root_norm , f"gamma{g_int}"),
            os.path.join("results", "morse", f"g{g_int}"),
            os.path.join("results", "morse", f"gamma{g_int}"),
        ]:
            candidates.append(os.path.join(d, "morse_potential.csv"))
    candidates += [
        os.path.join("results", "morse_potential.csv"),
        os.path.join("..", "results", "morse_potential.csv"),
        os.path.join("build", "Release", "morse_potential.csv"),
        os.path.join("build", "morse_potential.csv"),
    ]
    csv_path = next((c for c in candidates if os.path.exists(c)), None)
    if csv_path is None:
        for dirpath, _, filenames in os.walk(args.root):
            if "morse_potential.csv" in filenames:
                csv_path = os.path.join(dirpath, "morse_potential.csv"); break
    if csv_path is None:
        raise SystemExit("Nie znaleziono morse_potential.csv (sprawdź --root lub użyj --csv).")

    x, U = read_xy_csv(csv_path)

    # Poziomy energii
    Es = None
    if not args.no_levels and args.gamma is not None:
        nmax = morse_nmax(args.gamma)
        L = (nmax + 1) if args.levels is None else args.levels
        Es = morse_energies(args.gamma, L)

    # Zakresy
    (xl, xr), (yl, yr) = autodomain(
        x, U, xleft=args.xleft, xright=args.xright, xmax_sym=args.xmax,
        ymin=args.ymin, ymax=args.ymax, extra_y=Es
    )

    # Rysowanie
    plt.figure(figsize=(args.figw, args.figh), dpi=args.dpi)
    plt.plot(x, U, label="Morse U(x)")
    if args.grid: plt.grid(True, alpha=0.25)
    plt.axhline(0, lw=1, ls="--")


    if args.asymptote:
        y_as = 1.0
        plt.axhline(y_as, lw=1.2, ls="--", alpha=0.8)
        plt.text(
            xl + 0.02*(xr - xl), y_as - 0.02*(yr - yl),  # ЧУТЬ НИЖЕ линии
            f"asymptota = {y_as:g}",
            va="top", fontsize=9,
            bbox=dict(facecolor="white", alpha=0.7, pad=2, edgecolor="none")
        )


    # подпись количества уровней — переносим вправо от ямы
    if Es is not None and len(Es):
        x_pos = xl + 0.03 * (xr - xl)
        plt.text(
            x_pos, Es[-1],  # справа от «лесенки»
            f"poziomy: {len(Es)} (gamma={args.gamma:g})",
            va="bottom", fontsize=9,
            bbox=dict(facecolor="white", alpha=0.7, pad=2, edgecolor="none")
        )
        # Rysujemy poziomy energii jako przerywane linie, aby były widoczne na wykresie
        for e in Es:
            plt.hlines(e, xl, xr, colors="tab:orange", linestyles="--", linewidth=0.9, alpha=0.9)

    plt.xlim(xl, xr); plt.ylim(yl, yr)
    plt.title(args.title); plt.xlabel("x"); plt.ylabel("U(x)")
    plt.legend(loc="best"); plt.tight_layout()
    if args.save:
        os.makedirs(os.path.dirname(args.save), exist_ok=True)
        plt.savefig(args.save, dpi=max(200, args.dpi))
    plt.show()

if __name__ == "__main__":
    main()
