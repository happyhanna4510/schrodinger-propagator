#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Porównanie logów ewolucji:
- Oczekiwane kolumny: method,step,t,dt,dt_ms,matvecs,norm_err,theta_abs,theta_rel
- Wybór najlepszego przebiegu dla pary (gamma, dt) po minimalnym końcowym theta_abs.
- Dla Taylora wybieramy „najlepsze K” (wnioskowane z trybu matvecs) wśród dostępnych logów.
- Wykresy ograniczone wyłącznie do: norm_err oraz theta_abs (skala log dla błędów).
"""

import argparse, os, re, math
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext
import numpy as np


# --- Narzędzia IO ---

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def list_candidate_files(root: str):
    """Zbierz kandydackie pliki z logami (CSV/TSV/LOG)."""
    pats = ["**/*.csv","**/*.CSV","**/*.tsv","**/*.log","*.csv","*.CSV","*.tsv","*.log"]
    seen, out = set(), []
    for pat in pats:
        for f in glob(os.path.join(root, pat), recursive=True):
            f = os.path.abspath(f)
            if os.path.isfile(f) and f not in seen:
                seen.add(f); out.append(f)
    return sorted(out)

def sniff_sep(path: str) -> str:
    """Prosta heurystyka: wybierz najczęstszy z kandydatów jako separator."""
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        head = fh.read(4096)
    cands = [",",";","\t","|"]
    counts = {c: head.count(c) for c in cands}
    sep = max(counts, key=counts.get)
    return sep if counts[sep] > 0 else ","


# --- Ekstrakcja gamma z nazwy pliku/ścieżki ---

_GAMMA_RX = [
    re.compile(r"(?:^|[^\w])gamma[=_\-]?(\d+)(?:[^\w]|$)", re.I),
    re.compile(r"(?:^|[^\w])g(\d+)(?:[^\w]|$)", re.I),
    re.compile(r"_g(\d+)_", re.I),
]

def infer_gamma(df: pd.DataFrame, path: str) -> str:
    """1) jeśli kolumna gamma istnieje — użyj jej; 2) w przeciwnym razie spróbuj z nazwy pliku."""
    if "gamma" in df.columns:
        g = df["gamma"].mode(dropna=True)
        if len(g): return str(g.iloc[0])
    s = os.path.basename(path).replace("\\", "/")
    for rx in _GAMMA_RX:
        m = rx.search(s)
        if m: return m.group(1)
    m = re.search(r"g(\d+)(?=_[nN]\d+)", s)
    return m.group(1) if m else "unknown"


# --- Normalizacja metody ---

def canonical_method(name: str) -> str:
    n = (name or "").strip().lower()
    if n in ("taylor",): return "taylor"
    if n in ("cheb","chebyshev","chebyshev-polynomial","chebyshev_poly"): return "chebyshev"
    if n in ("rk4","rk-4","rk_4","runge-kutta4","rungekutta4"): return "rk4"
    return n  # jak jest


# --- Wczytanie i walidacja ---

NEEDED = ["method","step","t","dt","dt_ms","matvecs","norm_err","theta_abs","theta_rel"]

def load_one(path: str) -> pd.DataFrame:
    """Wczytaj jeden plik logu i przygotuj kolumny."""
    sep = sniff_sep(path)
    try:
        df = pd.read_csv(path, sep=sep, engine="python")
    except Exception as e:
        raise ValueError(f"read_csv nie powiodło się ({e})")
    missing = [c for c in NEEDED if c not in df.columns]
    if missing:
        raise ValueError(f"brak kolumn: {missing} (sep='{sep}')")

    # typy liczbowe (theta_rel zostaje tylko jako kolumna, ale nie rysujemy jej)
    num_cols = ["step","t","dt","dt_ms","matvecs","norm_err","theta_abs","theta_rel"]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # metoda, gamma, ścieżka
    df["method"] = df["method"].astype(str).map(canonical_method)
    df["gamma_str"] = infer_gamma(df, path)
    df["run_path"] = path

    # Szacunek K dla Taylora: tryb (mode) po matvecs; jeśli brak, użyj średniej
    if df["method"].iloc[0] == "taylor":
        try:
            k_mode = int(df["matvecs"].mode(dropna=True).iloc[0])
        except Exception:
            k_mode = int(round(df["matvecs"].dropna().mean())) if df["matvecs"].notna().any() else -1
        df["taylor_K"] = k_mode
    else:
        df["taylor_K"] = pd.NA

    return df


# --- Selekcja najlepszej trajektorii po minimalnym theta_abs na końcu ---

def _safe_pos_float(x) -> float:
    try:
        v = float(x)
        return v if (math.isfinite(v)) else math.inf
    except Exception:
        return math.inf

def pick_best_by_thetaabs_min(dfs):
    """
    Z kandydatów wybierz DF o najmniejszym końcowym theta_abs.
    Remis rozstrzygamy niższym średnim matvecs.
    """
    cand = []
    for dfi in dfs:
        if len(dfi) == 0: continue
        th_end = _safe_pos_float(dfi.iloc[-1].get("theta_abs", math.inf))
        if not math.isfinite(th_end): continue
        mv = float(dfi["matvecs"].dropna().mean()) if dfi["matvecs"].notna().any() else math.inf
        cand.append((th_end, mv, dfi))
    if not cand: return None
    cand.sort(key=lambda x: (x[0], x[1]))
    best_val = cand[0][0]
    near = [c for c in cand if c[0] <= 1.1 * best_val]
    near.sort(key=lambda x: (x[1], x[0]))
    return near[0][2]


# --- Rysowanie ---

ERR_METRICS = {"norm_err", "theta_abs"}  # tylko te dwie

def plot_metric(ax, df: pd.DataFrame, label: str, metric: str):
    d = df.sort_values("t").copy()
    if metric in ERR_METRICS:
        vals = d[metric].to_numpy()
        pos = vals[vals > 0]
        if pos.size == 0:
            return  # nic pozytywnego do log-skali
        vmin = pos.min()
        vals = np.where(vals > 0, vals, 0.1 * vmin)
        d[metric] = vals
        ax.plot(d["t"].values, d[metric].values, label=label)
        ax.set_yscale("log")
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_major_formatter(LogFormatterMathtext())
    else:
        ax.plot(d["t"].values, d[metric].values, label=label)
    ax.set_xlabel("t")
    ax.set_ylabel(metric)


# --- Główna logika ---

def main():
    ap = argparse.ArgumentParser(description="Porównanie logów metod (norm_err i theta_abs)")
    ap.add_argument("--root", required=True, help="Folder z logami (wyszukiwanie rekurencyjne).")
    ap.add_argument("--out", default="./plots_out", help="Folder wyjściowy na wykresy.")
    ap.add_argument("--dpi", type=int, default=140, help="DPI zapisanych obrazów.")
    ap.add_argument("--verbose", action="store_true", help="Więcej informacji diagnostycznych.")
    args = ap.parse_args()

    ensure_dir(args.out)
    files = list_candidate_files(args.root)
    if args.verbose:
        print(f"[INFO] Szukam plików w: {os.path.abspath(args.root)}")
        print(f"[INFO] Znaleziono kandydatów: {len(files)}")

    loaded, skipped = [], []
    for p in files:
        try:
            df = load_one(p)
            loaded.append(df)
        except Exception as e:
            skipped.append((p, str(e)))

    if args.verbose:
        print(f"[INFO] Poprawnie wczytane: {len(loaded)}; Pominięte: {len(skipped)}")
        for p, reason in skipped[:10]:
            print(f"[SKIP] {p}: {reason}")
        if len(skipped) > 10:
            print(f"[SKIP] ... jeszcze {len(skipped)-10} szt.")

    if not loaded:
        print("[WARN] Nie znaleziono żadnych poprawnych logów (sprawdź rozszerzenia i separator).")
        return

    big = pd.concat(loaded, ignore_index=True)

    # porządkujemy zestawy gamma i dt
    gammas = sorted(big["gamma_str"].unique(), key=lambda x: (x == "unknown", x))
    dts = sorted(big["dt"].dropna().unique())

    # kubełki (gamma, dt, method) -> lista DF (różne przebiegi dla tej samej pary)
    bucket = {}
    for df in loaded:
        if len(df) == 0: continue
        g = df["gamma_str"].iloc[0]
        dt_mode = df["dt"].mode(dropna=True)
        if len(dt_mode) == 0: continue
        dtv = float(dt_mode.iloc[0])
        m = df["method"].iloc[0]
        key = (g, dtv, m)
        bucket.setdefault(key, []).append(df)

    metrics = ["norm_err", "theta_abs"]  # tylko te dwie
    total_plots = 0

    for g in gammas:
        for dtv in dts:
            trio = {}

            # Taylor: wybór najlepszego K po minimalnym końcowym theta_abs
            dfs_t = bucket.get((g, dtv, "taylor"), [])
            best_t = pick_best_by_thetaabs_min(dfs_t)
            if best_t is not None:
                k_best = best_t["taylor_K"].iloc[0] if "taylor_K" in best_t.columns else None
                trio["taylor"] = (best_t, f"Taylor (best K={k_best})" if k_best is not None else "Taylor")

            # Chebyshev
            dfs_c = bucket.get((g, dtv, "chebyshev"), [])
            best_c = pick_best_by_thetaabs_min(dfs_c) if dfs_c else None
            if best_c is not None:
                trio["chebyshev"] = (best_c, "Chebyshev")

            # RK4
            dfs_r = bucket.get((g, dtv, "rk4"), [])
            best_r = pick_best_by_thetaabs_min(dfs_r) if dfs_r else None
            if best_r is not None:
                trio["rk4"] = (best_r, "RK4")

            if not trio:
                continue

            for metric in metrics:
                fig, ax = plt.subplots(figsize=(8, 4.5))
                empty = True
                for _, (dfp, label) in trio.items():
                    if metric in dfp.columns and dfp[metric].notna().any():
                        plot_metric(ax, dfp, label, metric)
                        empty = False
                if empty:
                    plt.close(fig)
                    continue

                ax.set_title(f"gamma={g}, dt={dtv:g} — {metric}")
                ax.legend(); ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fname = f"cmp_gamma{g}_dt{dtv:.0e}_{metric}.png".replace("+", "")
                fpath = os.path.join(args.out, fname)
                fig.savefig(fpath, dpi=args.dpi)
                plt.close(fig)
                total_plots += 1

            # Krótkie podsumowanie dla Taylora (jeśli jest)
            if "taylor" in trio:
                df_best = trio["taylor"][0]
                last = df_best.iloc[-1]
                th_final = _safe_pos_float(last.get("theta_abs"))
                print(f"[INFO] gamma={g} dt={dtv:g}: Taylor K={df_best['taylor_K'].iloc[0]} "
                      f"(final theta_abs={th_final:.3e}, mean matvecs={df_best['matvecs'].mean():.2f})")

    out_abs = os.path.abspath(args.out)
    if total_plots == 0:
        print(f"[WARN] Nie wygenerowano żadnych wykresów. Typowe powody:\n"
              f"  • nieprawidłowe nazwy metod (oczekiwane: cheb/chebyshev, rk4, taylor);\n"
              f"  • nietypowy separator — heurystyka zwykle działa, ale sprawdź kolumny;\n"
              f"  • brak dodatnich wartości w metrykach błędu (log-skala).")
    else:
        print(f"[OK] Zapisano wykresów: {total_plots} → {out_abs}")

    if skipped:
        print(f"[NOTE] Pominięte pliki: {len(skipped)} (użyj --verbose, aby zobaczyć szczegóły).")


if __name__ == "__main__":
    main()
