#!/usr/bin/env python3
"""
Plot Chebyshev expansion coefficients saved with --cheb-beta-log.
"""

import argparse
import csv
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt


def load_betas(csv_path: Path) -> Tuple[List[int], List[float]]:
    n_values: List[int] = []
    beta_abs: List[float] = []
    with csv_path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_values.append(int(row["n"]))
            beta_abs.append(float(row["beta_abs"]))
    return n_values, beta_abs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot |beta_n| vs n for Chebyshev expansion coefficients "
            "saved with --cheb-beta-log."
        )
    )
    parser.add_argument(
        "csv",
        type=Path,
        help="CSV file containing Chebyshev coefficients (columns: n, beta_abs, ...)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output image path (default: <csv-stem>_beta.png in the same directory)",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional plot title (default: derived from input filename)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.csv.exists():
        raise SystemExit(f"CSV file not found: {args.csv}")

    out_path = args.out
    if out_path is None:
        out_path = args.csv.with_name(f"{args.csv.stem}_beta.png")

    n_values, beta_abs = load_betas(args.csv)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(n_values, beta_abs, marker="o", linestyle="-", color="#1f77b4")
    ax.set_xlabel("n (Chebyshev order)")
    ax.set_ylabel(r"$|\beta_n|$")
    ax.set_title(args.title or f"Chebyshev coefficients from {args.csv.name}")
    ax.grid(True, which="both", linestyle=":", linewidth=0.8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")


if __name__ == "__main__":
    main()
