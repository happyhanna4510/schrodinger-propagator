import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_density(csv_path: Path):
    df = pd.read_csv(csv_path)
    if "t" not in df.columns:
        raise ValueError(f"column 't' not found in {csv_path}")
    x_cols = [c for c in df.columns if c != "t"]
    try:
        x = np.array([float(c) for c in x_cols], dtype=float)
    except ValueError as exc:
        raise ValueError(f"cannot parse x columns in {csv_path}") from exc
    t = df["t"].to_numpy(dtype=float)
    rho = df[x_cols].to_numpy(dtype=float)
    return t, x, rho


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def map_results_to_plots(results_dir: Path) -> Path:
    parts = results_dir.resolve().parts
    if "results" in parts:
        idx = parts.index("results")
        mapped = Path(*parts[:idx], "plots_out", *parts[idx + 1:])
    else:
        mapped = results_dir.parent / "plots_out" / results_dir.name
    return mapped


def plot_heatmap(data: np.ndarray, x: np.ndarray, t: np.ndarray, out_path: Path, title: str):
    plt.figure()
    plt.imshow(
        data,
        extent=[float(np.min(x)), float(np.max(x)), float(np.min(t)), float(np.max(t))],
        aspect="auto",
        origin="lower",
    )
    plt.xlabel("x")
    plt.ylabel("t")
    plt.title(title)
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Plot reference vs numerical densities. Reads num_density.csv and ref_density.csv "
            "from a simulation directory and writes heatmaps alongside diff."
        )
    )
    ap.add_argument(
        "--results-dir",
        required=True,
        help="Directory containing num_density.csv and ref_density.csv (e.g., results/U0_-1/g10/dt1e-5)",
    )
    ap.add_argument(
        "--out-dir",
        help="Optional output directory for PNGs. Default mirrors results_dir under plots_out.",
    )

    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    num_csv = results_dir / "num_density.csv"
    ref_csv = results_dir / "ref_density.csv"

    t_num, x_num, rho_num = load_density(num_csv)
    t_ref, x_ref, rho_ref = load_density(ref_csv)

    if not np.array_equal(x_num, x_ref):
        raise ValueError("x grids differ between numerical and reference density")
    if not np.array_equal(t_num, t_ref):
        raise ValueError("time grids differ between numerical and reference density")

    out_dir = Path(args.out_dir) if args.out_dir else map_results_to_plots(results_dir)
    ensure_dir(out_dir)

    heatmap_num = out_dir / "heatmap_num.png"
    heatmap_ref = out_dir / "heatmap_ref.png"
    heatmap_diff = out_dir / "heatmap_diff.png"

    plot_heatmap(rho_num, x_num, t_num, heatmap_num, "Numerical density |ψ|²")
    plot_heatmap(rho_ref, x_ref, t_ref, heatmap_ref, "Reference density |ψ_ref|²")
    plot_heatmap(rho_ref - rho_num, x_ref, t_ref, heatmap_diff, "ρ_ref - ρ_num")

    print(f"[OK] Saved {heatmap_num}")
    print(f"[OK] Saved {heatmap_ref}")
    print(f"[OK] Saved {heatmap_diff}")


if __name__ == "__main__":
    main()

