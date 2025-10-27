#!/usr/bin/env python3
"""
Universal comparator for time‑propagator logs (Taylor / RK4 / Chebyshev).

Assumptions
-----------
Each CSV has at least the columns:
  method, step, t, dt, dt_ms, matvecs, norm_err, theta
(Extra columns are ignored.)

File name convention (flexible):
  <method>[_K<k>|_m<m>|_p<p>]_dt<dt>_g<gamma>_N<...>_x<...>[...].csv
Examples:
  taylor_K6_dt1e-4_g10_N2001_x30.csv
  rk4_dt1e-5_g20_N2001_x30.csv
  chebyshev_m50_dt1e-6_g10_N2001_x30.csv

What it does
------------
- Scans a directory for CSVs.
- Parses parameters from filenames and augments the dataframe.
- Filters by (gamma, dt) so we compare apples to apples.
- Optionally auto-selects the "best" Taylor K (by lowest total runtime or matvecs).
- Plots, for the selected files, three figures:
    1) norm_err vs t
    2) theta vs t
    3) cumulative wall time vs t  (cum sum of dt_ms) and cumulative matvecs vs t
- Saves PNGs in --outdir.
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
from pathlib import Path

FNAME_RE = re.compile(
    r'(?P<method>[a-zA-Z0-9]+)'
    r'(?:_K(?P<K>\d+))?'
    r'(?:_m(?P<m>\d+))?'
    r'(?:_p(?P<p>\d+))?'
    r'_dt(?P<dt>[^_]+)'
    r'_g(?P<gamma>\d+)'
)

def parse_from_name(name: str):
    m = FNAME_RE.search(name)
    if not m:
        return {}
    d = m.groupdict()
    # normalize numeric fields
    for key in ('K','m','p','gamma'):
        if d.get(key):
            d[key] = int(d[key])
    # dt can be like '1e-4' or '0.0001'
    # keep as string for exact filtering, but also store numeric if possible
    try:
        d['dt_num'] = float(d['dt'].replace('e', 'e'))
    except Exception:
        d['dt_num'] = None
    return d

def load_one(path: str):
    base = os.path.basename(path)
    meta = parse_from_name(base)
    try:
        df = pd.read_csv(path)
    except Exception as e:
        raise RuntimeError(f"Failed reading {path}: {e}")
    # If columns missing in CSV, try to fill from filename
    for k,v in meta.items():
        if k not in df.columns:
            df[k] = v
    # Normalize required columns existence
    required = ['method','step','t','dt','dt_ms','matvecs','norm_err','theta']
    for col in required:
        if col not in df.columns:
            # Create if missing
            if col == 'method' and 'method' in meta:
                df['method'] = meta['method']
            elif col == 'dt' and 'dt' in meta:
                df['dt'] = meta['dt']
            elif col == 't':
                # try to infer from step * dt
                if 'step' in df.columns and 'dt' in df.columns:
                    try:
                        df['t'] = (df['step'] - df['step'].iloc[0]) * float(str(df['dt'].iloc[0]))
                    except Exception:
                        df['t'] = np.arange(len(df))
                else:
                    df['t'] = np.arange(len(df))
            else:
                # fallback zeros
                df[col] = 0
    df['label_extra'] = ''
    if 'K' in df.columns and pd.notna(df['K']).any():
        df['label_extra'] = 'K=' + df['K'].astype('Int64').astype(str).fillna('')
    elif 'm' in df.columns and pd.notna(df['m']).any():
        df['label_extra'] = 'm=' + df['m'].astype('Int64').astype(str).fillna('')
    elif 'p' in df.columns and pd.notna(df['p']).any():
        df['label_extra'] = 'p=' + df['p'].astype('Int64').astype(str).fillna('')
    df['fname'] = base
    return df, meta

def scan_dir(root: str, include_wide: bool=False):
    pats = ["*.csv"]
    files = []
    for pat in pats:
        files.extend(glob.glob(os.path.join(root, pat)))
    # filter out '_abs2_wide' unless requested
    if not include_wide:
        files = [f for f in files if "_abs2_wide" not in os.path.basename(f)]
    data = []
    metas = []
    for f in sorted(files):
        try:
            df, meta = load_one(f)
            data.append(df)
            metas.append(meta)
        except Exception as e:
            print(e)
    if not data:
        raise SystemExit("No CSV files found.")
    big = pd.concat(data, ignore_index=True)
    return big

def pick_best_taylor(df, metric="time"):
    """
    Given a subset dataframe (same gamma & dt), choose Taylor K that minimizes:
      - 'time' : total runtime (sum dt_ms)
      - 'matvecs' : total matvecs
      - 'final_err' : norm_err at max t
    """
    taylor = df[df['method'].str.lower() == 'taylor'].copy()
    if taylor.empty:
        return None

    def agg(g):
        out = {}
        out['tot_time_ms'] = g['dt_ms'].sum()
        out['tot_matvecs'] = g['matvecs'].sum()
        # final error at largest t in this group
        idx = g['t'].idxmax()
        out['final_err'] = g.loc[idx, 'norm_err']
        out['K'] = g['K'].dropna().iloc[0] if 'K' in g else None
        out['fname'] = g['fname'].iloc[0]
        return pd.Series(out)

    tab = taylor.groupby('fname', group_keys=False).apply(agg).reset_index(drop=True)

    if metric == "time":
        row = tab.sort_values('tot_time_ms', ascending=True).iloc[0]
    elif metric == "matvecs":
        row = tab.sort_values('tot_matvecs', ascending=True).iloc[0]
    elif metric == "final_err":
        row = tab.sort_values('final_err', ascending=True).iloc[0]
    else:
        row = tab.iloc[0]

    chosen_fname = row['fname']
    best = taylor[taylor['fname'] == chosen_fname]
    return best

def plot_lines(ax, df, label):
    ax.plot(df['t'].values, df['y'].values, label=label)
    ax.set_xlabel("t")
    ax.legend()
    ax.grid(True, which="both", linestyle=":")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Directory with CSV logs")
    ap.add_argument("--gamma", required=True, type=int, help="Select g (e.g., 10 or 20)")
    ap.add_argument("--dt", required=True, help="Select dt string exactly as in filename, e.g., 1e-4, 1e-5, 1e-6")
    ap.add_argument("--outdir", default="plots_out", help="Where to save PNGs")
    ap.add_argument("--include-wide", action="store_true", help="Include *_abs2_wide.csv files")
    ap.add_argument("--taylor-pick", choices=["time","matvecs","final_err","all"], default="time",
                    help="How to pick single Taylor K for overlay; 'all' keeps all K")
    ap.add_argument("--also_save_csv_summary", action="store_true", help="Save a per-file summary CSV")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    big = scan_dir(args.root, include_wide=args.include_wide)

    # Keep only selected dt & gamma
    # gamma may be in column 'gamma' or parsed via name (already placed in df['gamma'])
    sel = big[
        (np.isclose(big['gamma'].astype(float), float(args.gamma), rtol=1e-8))
        & (np.isclose(big['dt'].astype(float), float(args.dt), rtol=1e-8))
    ].copy()
    if sel.empty:
        raise SystemExit(f"No rows match gamma={args.gamma} and dt={args.dt}.")

    # Optional summary
    if args.also_save_csv_summary:
        def summarize(g):
            return pd.Series({
                "method": g['method'].iloc[0],
                "K": g['K'].dropna().iloc[0] if 'K' in g and not g['K'].dropna().empty else None,
                "tot_time_ms": g['dt_ms'].sum(),
                "tot_matvecs": g['matvecs'].sum(),
                "final_err": g.loc[g['t'].idxmax(), 'norm_err'],
                "rows": len(g),
                "fname": g['fname'].iloc[0],
            })
        summ = sel.groupby('fname').apply(summarize).reset_index(drop=True)
        summ_out = os.path.join(args.outdir, f"summary_g{args.gamma}_dt{args.dt}.csv")
        summ.to_csv(summ_out, index=False)
        print(f"Saved summary: {summ_out}")

    # Optionally collapse Taylor to one K
    parts = []
    if args.taylor_pick == "all":
        parts.append(sel[sel['method'].str.lower()=='taylor'])
    else:
        best_taylor = pick_best_taylor(sel, metric=args.taylor_pick)
        if best_taylor is not None:
            parts.append(best_taylor)
    # Non-Taylor methods
    parts.append(sel[sel['method'].str.lower()!='taylor'])
    data = pd.concat(parts, ignore_index=True) if parts else sel

    # Build plotting groups (each file = one curve)
    groups = list(data.groupby('fname'))

    # 1) norm_err vs t
    fig1, ax1 = plt.subplots(figsize=(8,5))
    for fname, g in groups:
        label = f"{g['method'].iloc[0]}"
        extra = g['label_extra'].dropna().iloc[0] if 'label_extra' in g and not g['label_extra'].dropna().empty else ''
        if extra:
            label += f" ({extra})"
        dfp = g[['t','norm_err']].rename(columns={'norm_err':'y'}).sort_values('t')
        plot_lines(ax1, dfp, label=label)
    ax1.set_ylabel("norm_err")
    ax1.set_yscale("log")
    fig1.tight_layout()
    out1 = os.path.join(args.outdir, f"norm_err_g{args.gamma}_dt{args.dt}.png")
    fig1.savefig(out1, dpi=150)

    # 2) theta vs t
    fig2, ax2 = plt.subplots(figsize=(8,5))
    for fname, g in groups:
        label = f"{g['method'].iloc[0]}"
        extra = g['label_extra'].dropna().iloc[0] if 'label_extra' in g and not g['label_extra'].dropna().empty else ''
        if extra:
            label += f" ({extra})"
        dfp = g[['t','theta']].rename(columns={'theta':'y'}).sort_values('t')
        plot_lines(ax2, dfp, label=label)
    ax2.set_ylabel("theta")
    ax2.set_yscale("log")
    fig2.tight_layout()
    out2 = os.path.join(args.outdir, f"theta_g{args.gamma}_dt{args.dt}.png")
    fig2.savefig(out2, dpi=150)

    # 3) cumulative time & matvecs vs t
    fig3, ax3 = plt.subplots(figsize=(8,5))
    for fname, g in groups:
        label = f"{g['method'].iloc[0]}"
        extra = g['label_extra'].dropna().iloc[0] if 'label_extra' in g and not g['label_extra'].dropna().empty else ''
        if extra:
            label += f" ({extra})"
        gg = g.sort_values('t').copy()
        gg['cum_ms'] = gg['dt_ms'].cumsum()
        dfp = gg[['t','cum_ms']].rename(columns={'cum_ms':'y'})
        plot_lines(ax3, dfp, label=label)
    ax3.set_ylabel("cumulative time [ms]")
    fig3.tight_layout()
    out3 = os.path.join(args.outdir, f"cumtime_g{args.gamma}_dt{args.dt}.png")
    fig3.savefig(out3, dpi=150)

    fig4, ax4 = plt.subplots(figsize=(8,5))
    for fname, g in groups:
        label = f"{g['method'].iloc[0]}"
        extra = g['label_extra'].dropna().iloc[0] if 'label_extra' in g and not g['label_extra'].dropna().empty else ''
        if extra:
            label += f" ({extra})"
        gg = g.sort_values('t').copy()
        gg['cum_mv'] = gg['matvecs'].cumsum()
        dfp = gg[['t','cum_mv']].rename(columns={'cum_mv':'y'})
        plot_lines(ax4, dfp, label=label)
    ax4.set_ylabel("cumulative matvecs")
    fig4.tight_layout()
    out4 = os.path.join(args.outdir, f"cummatvecs_g{args.gamma}_dt{args.dt}.png")
    fig4.savefig(out4, dpi=150)

    print(f"Saved plots to: {args.outdir}")
    print(out1)
    print(out2)
    print(out3)
    print(out4)

if __name__ == "__main__":
    main()
