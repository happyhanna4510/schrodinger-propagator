# scripts/plot_psi_wide_all.py
import argparse, os, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def load_wide(path: str):
    df = pd.read_csv(path, dtype=str, low_memory=False)
    t = pd.to_numeric(df['t'], errors='coerce').astype(float).to_numpy()
    x_cols = df.columns[1:]
    x = np.array([float(c) for c in x_cols])
    Y = df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').to_numpy()
    return t, x, Y

def pick_indices(n, count=5):
    if n <= count: return list(range(n))
    return sorted({0, n//4, n//2, 3*n//4, n-1})

def plot_snapshots(path_in: str, title: str, ylabel: str, outpath: str, count=5):
    t, x, Y = load_wide(path_in)
    idxs = pick_indices(len(t), count=count)
    plt.figure()
    for i in idxs:
        plt.plot(x, Y[i], label=f"t={t[i]:.2e}")
    plt.xlabel("x"); plt.ylabel(ylabel); plt.title(title); plt.legend()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=160, bbox_inches="tight")
    plt.close()
    print(f"[saved] {outpath}")

def process_one(stem: str, results_dir: str, figures_dir: str, count: int):
    abs2 = os.path.join(results_dir, f"{stem}_abs2_wide.csv")
    re   = os.path.join(results_dir, f"{stem}_re_wide.csv")
    im   = os.path.join(results_dir, f"{stem}_im_wide.csv")

    if os.path.exists(abs2):
        plot_snapshots(abs2, f"|ψ(x,t)|² — {stem}", r"|ψ(x,t)|²",
                       os.path.join(figures_dir, f"{stem}_abs2_snapshots.png"),
                       count=count)
    else:
        print(f"[warn] not found: {abs2}")

    if os.path.exists(re):
        plot_snapshots(re, f"Re ψ(x,t) — {stem}", "Re ψ",
                       os.path.join(figures_dir, f"{stem}_re_snapshots.png"),
                       count=count)

    if os.path.exists(im):
        plot_snapshots(im, f"Im ψ(x,t) — {stem}", "Im ψ",
                       os.path.join(figures_dir, f"{stem}_im_snapshots.png"),
                       count=count)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results", help="папка с wide-CSV")
    ap.add_argument("--figures", default="figures", help="куда сохранить PNG")
    ap.add_argument("--count", type=int, default=5, help="сколько срезов по времени рисовать на графике")
    args = ap.parse_args()

    # Находим все abs2_wide и из них получаем stem
    pattern = os.path.join(args.results, "*_abs2_wide.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[warn] no files: {pattern}")
        return

    for f in files:
        base = os.path.basename(f)
        stem = base.replace("_abs2_wide.csv", "")
        process_one(stem, args.results, args.figures, args.count)

    print("[done]")

if __name__ == "__main__":
    main()
