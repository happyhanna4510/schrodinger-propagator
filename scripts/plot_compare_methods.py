#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re, math
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext
import numpy as np


def ensure_dir(p): os.makedirs(p, exist_ok=True)

def list_candidate_files(root: str):
    pats = ["**/*.csv","**/*.CSV","**/*.tsv","**/*.log","*.csv","*.CSV","*.tsv","*.log"]
    seen = set()
    out = []
    for pat in pats:
        for f in glob(os.path.join(root, pat), recursive=True):
            f = os.path.abspath(f)
            if os.path.isfile(f) and f not in seen:
                seen.add(f); out.append(f)
    return sorted(out)

# авторазделитель
def sniff_sep(path: str):
    import io
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        head = fh.read(4096)
    # простая эвристика: чаще встречающийся из кандидатов
    cands = [",",";","\t","|"]
    counts = {c: head.count(c) for c in cands}
    sep = max(counts, key=counts.get)
    # если вообще нет разделителей — пусть будет запятая (потом отфильтруем по колонкам)
    return sep if counts[sep] > 0 else ","

_gamma_rx = [
    re.compile(r"gamma[=_\-]?(\d+)", re.I),
    re.compile(r"\bg(\d+)\b", re.I),
]

# --- gamma from path/file name ---
_gamma_rx = [
    re.compile(r"(?:^|[^\w])gamma[=_\-]?(\d+)(?:[^\w]|$)", re.I),
    re.compile(r"(?:^|[^\w])g(\d+)(?:[^\w]|$)", re.I),          # ловит _g10_, -g20, g10.
    re.compile(r"_g(\d+)_", re.I),                               # узкий случай: _g10_
]

def infer_gamma(df: pd.DataFrame, path: str):
    # 1) если есть колонка — используем её
    if "gamma" in df.columns:
        g = df["gamma"].mode(dropna=True)
        if len(g): return str(g.iloc[0])
    # 2) иначе — берём из имени файла/пути
    s = os.path.basename(path).replace("\\", "/")
    for rx in _gamma_rx:
        m = rx.search(s)
        if m:
            return m.group(1)
    # 3) как крайний случай — ищем 'gNN' прямо перед 'N2001' и т.п.
    m = re.search(r"g(\d+)(?=_[nN]\d+)", s)
    return m.group(1) if m else "unknown"


# приведение имён методов к канону
def canonical_method(name: str):
    n = (name or "").strip().lower()
    if n in ("taylor",): return "taylor"
    if n in ("cheb","chebyshev","chebyshev-polynomial","chebyshev_poly"): return "chebyshev"
    if n in ("rk4","rk-4","rk_4","runge-kutta4","rungekutta4"): return "rk4"
    return n  # как есть

NEEDED = ["method","step","t","dt","dt_ms","matvecs","norm_err","theta","e_true"]

def safe_float(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else math.inf
    except Exception:
        return math.inf

def load_one(path: str):
    sep = sniff_sep(path)
    try:
        df = pd.read_csv(path, sep=sep, engine="python")
    except Exception as e:
        raise ValueError(f"read_csv failed ({e})")
    missing = [c for c in NEEDED if c not in df.columns]
    if missing: raise ValueError(f"missing cols: {missing} (sep='{sep}')")
    # типы
    num_cols = ["step","t","dt","dt_ms","matvecs","norm_err","theta","e_true"]
    for c in num_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["method"] = df["method"].astype(str).map(canonical_method)
    df["gamma_str"] = infer_gamma(df, path)
    df["run_path"] = path
    # оценка K по mode(matvecs) для тэйлора
    if df["method"].iloc[0] == "taylor":
        try:
            k_mode = int(df["matvecs"].mode(dropna=True).iloc[0])
        except Exception:
            k_mode = int(round(df["matvecs"].dropna().mean())) if df["matvecs"].notna().any() else -1
        df["taylor_K"] = k_mode
    else:
        df["taylor_K"] = pd.NA
    return df

def pick_best_by_etruemin(dfs):
    cand = []
    for dfi in dfs:
        if len(dfi)==0: continue
        e_end = abs(safe_float(dfi.iloc[-1].get("e_true", math.inf)))
        if not math.isfinite(e_end): continue
        mv = float(dfi["matvecs"].dropna().mean()) if dfi["matvecs"].notna().any() else math.inf
        cand.append((e_end, mv, dfi))
    if not cand: return None
    cand.sort(key=lambda x:(x[0],x[1]))
    best_e = cand[0][0]
    near = [c for c in cand if c[0] <= 1.1*best_e]
    near.sort(key=lambda x:(x[1],x[0]))
    return near[0][2]

def plot_metric(ax, df, label, metric):
    d = df.sort_values("t").copy()

    # Логарифмические метрики — рисуем в log10
    ERR_METRICS = {"norm_err", "theta", "e_true"}

    if metric in ERR_METRICS:
        # оставим только положительные; нули заменим на минимально-положительное
        vals = d[metric].to_numpy()
        # вычислим минимально положительное
        pos = vals[vals > 0]
        if pos.size == 0:
            return  # нечего рисовать
        vmin = pos.min()
        # заменим неположительные на 0.1*vmin для корректного лог-масштаба
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="./plots_out")
    ap.add_argument("--dpi", type=int, default=140)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    ensure_dir(args.out)
    files = list_candidate_files(args.root)
    if args.verbose:
        print(f"[INFO] Поиск файлов в: {os.path.abspath(args.root)}")
        print(f"[INFO] Найдено кандидатов: {len(files)}")

    loaded, skipped = [], []
    for p in files:
        try:
            df = load_one(p)
            loaded.append(df)
        except Exception as e:
            skipped.append((p, str(e)))
    if args.verbose:
        print(f"[INFO] Прочитано корректных файлов: {len(loaded)}; Отброшено: {len(skipped)}")
        for p,reason in skipped[:10]:
            print(f"[SKIP] {p}: {reason}")
        if len(skipped)>10:
            print(f"[SKIP] ... ещё {len(skipped)-10} шт.")

    if not loaded:
        print("[WARN] Не найдено ни одного корректного лога (проверь расширения и разделители).")
        return

    big = pd.concat(loaded, ignore_index=True)
    gammas = sorted(big["gamma_str"].unique(), key=lambda x:(x=="unknown", x))
    dts = sorted(big["dt"].dropna().unique())
    metrics = ["dt_ms","matvecs","norm_err","theta","e_true"]

    # сгруппируем по файлам (gamma, dt, method)
    bucket = {}
    for df in loaded:
        if len(df)==0: continue
        g = df["gamma_str"].iloc[0]
        dt_mode = df["dt"].mode(dropna=True)
        if len(dt_mode)==0: continue
        dtv = float(dt_mode.iloc[0])
        m = df["method"].iloc[0]
        key = (g, dtv, m)
        bucket.setdefault(key, []).append(df)

    total_plots = 0
    for g in gammas:
        for dtv in dts:
            trio = {}
            # Taylor
            dfs_t = bucket.get((g, dtv, "taylor"), [])
            best_t = pick_best_by_etruemin(dfs_t)
            if best_t is not None:
                k_best = best_t["taylor_K"].iloc[0] if "taylor_K" in best_t.columns else None
                trio["taylor"] = (best_t, f"Taylor (best K={k_best})")
            # Chebyshev
            dfs_c = bucket.get((g, dtv, "chebyshev"), [])
            if dfs_c:
                best_c = pick_best_by_etruemin(dfs_c)
                trio["chebyshev"] = (best_c, "Chebyshev")
            # RK4
            dfs_r = bucket.get((g, dtv, "rk4"), [])
            if dfs_r:
                best_r = pick_best_by_etruemin(dfs_r)
                trio["rk4"] = (best_r, "RK4")

            if not trio: 
                continue

            for metric in metrics:
                fig, ax = plt.subplots(figsize=(8,4.5))
                empty = True
                for _, (dfp, label) in trio.items():
                    if metric in dfp.columns and dfp[metric].notna().any():
                        plot_metric(ax, dfp, label, metric)
                        empty = False
                if empty:
                    plt.close(fig); continue
                ax.set_title(f"gamma={g}, dt={dtv:g} — {metric}")
                ax.legend(); ax.grid(True, alpha=0.3)
                fig.tight_layout()
                fname = f"cmp_gamma{g}_dt{dtv:.0e}_{metric}.png".replace("+","")
                fpath = os.path.join(args.out, fname)
                fig.savefig(fpath, dpi=args.dpi)
                plt.close(fig)
                total_plots += 1

            # консольная сводка
            if "taylor" in trio:
                df_best = trio["taylor"][0]
                last = df_best.iloc[-1]
                print(f"[INFO] gamma={g} dt={dtv:g}: Taylor K={df_best['taylor_K'].iloc[0]} "
                      f"(final |e_true|={abs(safe_float(last['e_true'])):.3e}, "
                      f"mean matvecs={df_best['matvecs'].mean():.2f})")

    out_abs = os.path.abspath(args.out)
    if total_plots==0:
        print(f"[WARN] Графики не построены. Частые причины:\n"
              f"  • метод назван не как 'cheb/chebyshev' или 'rk4/rk-4/rk_4';\n"
              f"  • нестандартный разделитель — теперь он авто-определяется, но проверь колонки;\n"
              f"  • логи для пар (gamma, dt) есть только у одного метода — это ок, но проверь, что колонки числовые.\n"
              f"Попробуй запустить с флагом --verbose для диагностики.")
    else:
        print(f"[OK] Сохранено графиков: {total_plots} → {out_abs}")
    if skipped:
        print(f"[NOTE] Отброшено файлов: {len(skipped)} (см. --verbose для деталей).")

if __name__ == "__main__":
    main()
