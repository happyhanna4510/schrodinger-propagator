import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

OUT_DIR_NAME = "plots_out"

def load_abs2_wide(csv_path: str):
    df = pd.read_csv(csv_path)
    if "t" not in df.columns:
        raise ValueError(f"'t' column not found in {csv_path}")
    x_cols = [c for c in df.columns if c != "t"]
    try:
        x = np.array([float(c) for c in x_cols], dtype=float)
    except ValueError as e:
        raise ValueError(f"Cannot parse x columns as floats in {csv_path}: {e}")
    t = df["t"].to_numpy(dtype=float)
    psi2 = df[x_cols].to_numpy(dtype=float)
    return t, x, psi2


def plot_heatmap(t, x, psi2, out_png):
    plt.figure()
    plt.imshow(
        psi2,
        extent=[np.min(x), np.max(x), np.min(t), np.max(t)],
        aspect="auto",
        origin="lower",
    )
    plt.xlabel("x")
    plt.ylabel("t")
    plt.title("|ψ(x,t)|² heatmap")
    plt.colorbar(label="|ψ|²")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_slices(t, x, psi2, out_png, nslices=5):
    plt.figure()
    if len(t) == 0:
        plt.close()
        return
    idx = np.unique(np.linspace(0, len(t) - 1, num=max(1, nslices)).astype(int))
    for k in idx:
        plt.plot(x, psi2[k, :], label=f"t={t[k]:.3f}")
    plt.xlabel("x")
    plt.ylabel("|ψ|²")
    plt.title("|ψ|² slices over time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def out_paths_for(csv_path: str, out_dir: str | None):
    stem = os.path.splitext(os.path.basename(csv_path))[0]
    if out_dir:
        ensure_dir(out_dir)
        heatmap_png = os.path.join(out_dir, f"{stem}__heatmap.png")
        slices_png  = os.path.join(out_dir, f"{stem}__slices.png")
    else:
        heatmap_png = f"{os.path.splitext(csv_path)[0]}__heatmap.png"
        slices_png  = f"{os.path.splitext(csv_path)[0]}__slices.png"
    return heatmap_png, slices_png


def process_file(csv_path: str, nslices: int, out_dir: str | None, dry_run: bool=False):
    try:
        t, x, psi2 = load_abs2_wide(csv_path)
    except Exception as e:
        print(f"[ERROR] {csv_path}: {e}", file=sys.stderr)
        return

    heatmap_png, slices_png = out_paths_for(csv_path, out_dir)

    if dry_run:
        print(f"[DRY] Read {csv_path}: t={len(t)} x={len(x)} array={psi2.shape}")
        print(f"[DRY] Save -> {heatmap_png}, {slices_png}")
        return

    plot_heatmap(t, x, psi2, heatmap_png)
    plot_slices(t, x, psi2, slices_png, nslices=nslices)
    print(f"[OK] {csv_path} -> {heatmap_png}, {slices_png}")


def discover_files(inputs, pattern, results_dir: str | None, script_dir: str):
    """
    Discovery rules tailored for project layout (lowercase):
    - If inputs provided -> use them verbatim.
    - Else if results_dir provided -> search there by pattern.
    - Else try CWD by pattern.
    - Else try CWD/results.
    - Else try <script_dir>/../results (common when script lives in Scripts/).
    - Else try <script_dir> (if someone dropped CSVs next to the script).
    """
    if inputs:
        return inputs

    candidates = []
    if results_dir:
        candidates.append(os.path.abspath(results_dir))

    cwd = os.getcwd()
    candidates.append(cwd)
    candidates.append(os.path.join(cwd, "results"))
    candidates.append(os.path.abspath(os.path.join(script_dir, os.pardir, "results")))
    candidates.append(script_dir)

    for base in candidates:
        if not base or not os.path.isdir(base):
            continue
        files = sorted(glob.glob(os.path.join(base, pattern)))
        if files:
            return files

    return []


def choose_default_out_dir(results_dir: str | None, script_dir: str) -> str:
    """
    Default output directory when --out-dir is not provided.
    Priority (lowercase folders):
      1) If results_dir is given -> <parent_of_results_dir>/plots_out
      2) If ./results exists in CWD -> ./plots_out
      3) If ../results exists relative to script -> ../plots_out
      4) Else -> ./plots_out (in CWD)
    """
    cwd = os.getcwd()

    # 1) Явно задан results_dir -> кладём в соседа 'plots_out'
    if results_dir:
        results_abs = os.path.abspath(results_dir)
        parent = os.path.dirname(results_abs)
        return os.path.join(parent, OUT_DIR_NAME)

    # 2) В текущем каталоге есть results -> кладём в ./plots_out
    if os.path.isdir(os.path.join(cwd, "results")):
        return os.path.join(cwd, OUT_DIR_NAME)

    # 3) Если скрипт лежит в Scripts/, а рядом есть ../results -> кладём в ../plots_out
    repo_results = os.path.abspath(os.path.join(script_dir, os.pardir, "results"))
    if os.path.isdir(repo_results):
        parent = os.path.dirname(repo_results)
        return os.path.join(parent, OUT_DIR_NAME)

    # 4) По умолчанию — ./plots_out
    return os.path.join(cwd, OUT_DIR_NAME)



def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    ap = argparse.ArgumentParser(
        description=(
            "Plot heatmap and slices from *_abs2_wide.csv files (|ψ|²).\n"
            "Designed for repo layout: Scripts/plot_abs2_wide.py and results/*.csv.\n"
            "Run from project root or anywhere; use --results-dir/--out-dir if needed.\n"
            "Default output dir: sibling 'plots_out' next to 'results' (or ./plots_out)."
        )
    )
    ap.add_argument(
        "inputs",
        nargs="*",
        help="Explicit CSV files to process. If omitted, discovery rules/pattern are used."
    )
    ap.add_argument(
        "--pattern",
        default="*_abs2_wide.csv",
        help="Glob pattern for discovery (default: *_abs2_wide.csv)."
    )
    ap.add_argument(
        "--results-dir",
        default=None,
        help="Directory to search for CSVs (e.g., results). Overrides discovery."
    )
    ap.add_argument(
        "--out-dir",
        default=None,
        help="Directory to save plots. Default: sibling 'plots_out' next to 'results' (or ./plots_out)."
    )
    ap.add_argument(
        "--slices",
        type=int,
        default=5,
        help="How many time slices to overlay (default: 5)."
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without creating plots."
    )
    args = ap.parse_args()

    files = discover_files(args.inputs, args.pattern, args.results_dir, script_dir)

    if not files:
        print(
            "No input files found.\n"
            "Try: python Scripts/plot_abs2_wide.py --results-dir results\n"
            ' or: python Scripts/plot_abs2_wide.py --pattern "*_abs2_wide.csv"',
            file=sys.stderr
        )
        sys.exit(1)

    # <-- ВАЖНО: вычисляем дефолтный каталог вывода
    out_dir = args.out_dir or choose_default_out_dir(args.results_dir, script_dir)
    ensure_dir(out_dir)

    for fp in files:
        process_file(fp, nslices=args.slices, out_dir=out_dir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
