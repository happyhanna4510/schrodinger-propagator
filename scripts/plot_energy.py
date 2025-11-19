"""
Plot energy E(t) from an energy CSV file produced with --log-energy.
This script never overwrites existing plots: if the target path exists,
"_N" suffixes are added automatically.
"""

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd


def make_unique_path(path: Path) -> Path:
    candidate = path
    counter = 1
    while candidate.exists():
        candidate = path.with_name(f"{path.stem}_{counter}{path.suffix}")
        counter += 1
    return candidate


def derive_output_path(csv_path: Path, out: Optional[str]) -> Path:
    target_dir = Path(out) if out else csv_path.parent
    target_dir.mkdir(parents=True, exist_ok=True)
    base = csv_path.stem + "_energy.png"
    candidate = target_dir / base
    return make_unique_path(candidate)


def plot_energy(csv_path: Path, out_path: Path, dpi: int) -> Path:
    df = pd.read_csv(csv_path)
    if not {"t", "E"}.issubset(df.columns):
        raise ValueError("CSV must contain 't' and 'E' columns")

    fig, ax = plt.subplots()
    ax.plot(df["t"], df["E"], label="E(t)")
    ax.set_xlabel("t")
    ax.set_ylabel("E")
    ax.set_title("Energy vs time")
    ax.grid(True, alpha=0.3)
    ax.legend()

    out_path = make_unique_path(out_path)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot energy vs time from CSV")
    ap.add_argument("csv", type=Path, help="Energy CSV with columns t,E")
    ap.add_argument("--out", type=str, help="Output directory or filename for the plot")
    ap.add_argument("--dpi", type=int, default=140, help="DPI of the saved plot")
    args = ap.parse_args()

    csv_path = args.csv
    if not csv_path.is_file():
        raise SystemExit(f"Input CSV not found: {csv_path}")

    out_path = Path(args.out) if args.out else derive_output_path(csv_path, None)
    if out_path.suffix.lower() not in {".png", ".pdf", ".svg"}:
        out_path = derive_output_path(csv_path, args.out)

    final_path = plot_energy(csv_path, out_path, args.dpi)
    print(f"Saved energy plot to: {final_path}")


if __name__ == "__main__":
    main()
