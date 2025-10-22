import argparse
import glob
import os
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def extract_K_from_name(path: str):
    """Return K as string (e.g., '4') or None if not present."""
    base = os.path.basename(path)
    mK = re.search(r'(?:^|[_\-])K(\d+)(?:[_\-]|\.|$)', base)
    return mK.group(1) if mK else None


def extract_dt_from_name(path: str):
    base = os.path.basename(path)
    mdt = (re.search(r'dt([0-9eE\-\.+]+)', base)
           or re.search(r'([0-9]+e-[0-9]+)', base))
    return mdt.group(1) if mdt else None


def guess_label(path: str) -> str:
    """Legend label: prefer K and dt."""
    K  = extract_K_from_name(path)
    dt = extract_dt_from_name(path)
    if K and dt: return f"K={K}, dt={dt}"
    if K:        return f"K={K}"
    if dt:       return f"dt={dt}"
    return os.path.basename(path)


def load_frames(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p)
        req = {'t','norm','p0','err_exact'}
        if not req.issubset(df.columns):
            raise ValueError(f"{p} must have columns: {sorted(req)}")
        df['_label'] = guess_label(p)
        df['_K']     = extract_K_from_name(p) or "unknown"
        df['_dt']    = extract_dt_from_name(p) or "unknown"
        df['_src']   = p
        frames.append(df)
    return frames



def _safe(s: str) -> str:
    return str(s).replace('+','').replace('-','m').replace('.','p')

def save_plots_by_dt(frames, outdir: Path, show: bool):
    """
    Делает графики err_exact(t) по dt: на каждом рисунке все K.
    Сохраняет: err_exact_dt_<dt>.png
    """
    groups = {}
    for df in frames:
        # взять СКАЛЯР, а не серию
        dt = df['_dt'].iloc[0] if '_dt' in df.columns else None
        if not dt:
            base = os.path.basename(df['_src']) if '_src' in df.columns else ''
            mdt = (re.search(r'dt([0-9eE\-\.+]+)', base)
                   or re.search(r'([0-9]+e-[0-9]+)', base))
            dt = mdt.group(1) if mdt else "unknown"
        dt = str(dt)
        groups.setdefault(dt, []).append(df)

    for dt, dfs in groups.items():
        plt.figure()
        for df in dfs:
            k = df['_K'].iloc[0] if '_K' in df.columns else "?"
            lab = f"K={k}"
            plt.plot(df['t'].values, df['err_exact'].values, label=lab)
        plt.xlabel("t")
        plt.ylabel("err_exact")
        plt.yscale("log")
        plt.legend()
        plt.title(f"err_exact(t) — dt={dt}")
        plt.savefig(outdir / f"err_exact_dt_{_safe(dt)}.png", dpi=160)
        if show:
            plt.show()
        else:
            plt.close()


def save_plots_for_group(frames, outdir: Path, k_tag: str, show: bool):
    """
    frames: list of dataframes (all with the same _K)
    outdir: where to save
    k_tag : e.g. 'K4' or 'Kunknown'
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # 1) norm(t)
    plt.figure()
    for df in frames:
        plt.plot(df['t'].values, df['norm'].values, label=df['_label'].iloc[0])
    plt.xlabel("t")
    plt.ylabel("norm")
    plt.legend()
    plt.title(f"norm(t) — {k_tag}")
    plt.savefig(outdir / f"norm_{k_tag}.png", dpi=160)

    # 2) p0(t)
    plt.figure()
    for df in frames:
        plt.plot(df['t'].values, df['p0'].values, label=df['_label'].iloc[0])
    plt.xlabel("t")
    plt.ylabel("p0")
    plt.legend()
    plt.title(f"p0(t) — {k_tag}")
    plt.savefig(outdir / f"p0_{k_tag}.png", dpi=160)

    # 3) err_exact(t)
    plt.figure()
    for df in frames:
        plt.plot(df['t'].values, df['err_exact'].values, label=df['_label'].iloc[0])
    plt.xlabel("t")
    plt.ylabel("err_exact")
    plt.yscale("log")
    plt.legend()
    plt.title(f"err_exact(t) — {k_tag}")
    plt.savefig(outdir / f"err_exact_{k_tag}.png", dpi=160)

    if show:
        plt.show()
    else:
        plt.close('all')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", action="append", default=[], help="Path to a taylor_log CSV (repeatable)")
    ap.add_argument("--glob", default="", help="Glob to collect CSVs, e.g. 'results/taylor_log*.csv'")
    ap.add_argument("--outdir", default="", help="Output directory (default = dir of the first CSV)")
    ap.add_argument("--no_show", action="store_true", help="Do not show windows; only save PNGs")
    ap.add_argument("--split_by_K", action="store_true",
                    help="If set, create subfolders K*/ even when only one K is present")
    args = ap.parse_args()

    paths = list(args.csv)
    if args.glob:
        paths.extend(sorted(glob.glob(args.glob)))
    if not paths:
        ap.error("No CSVs provided. Use --csv ... or --glob 'pattern'")

    frames = load_frames(paths)

    # базовый outdir
    outdir = Path(args.outdir) if args.outdir else Path("figures")
    outdir.mkdir(parents=True, exist_ok=True)


    # сгруппировать по K
    groups = {}
    for df in frames:
        groups.setdefault(df['_K'].iloc[0], []).append(df)

    # если K всего один и split_by_K не включён — сохраняем в базовую папку с суффиксом _K*
    if len(groups) == 1 and not args.split_by_K:
        only_k = next(iter(groups))
        k_tag = f"K{only_k}"
        save_plots_for_group(groups[only_k], outdir, k_tag, show=not args.no_show)
    else:
        # несколько K → делаем подпапки K4/, K6/, ...
        for k, dfs in groups.items():
            k_tag = f"K{k}"
            sub = outdir / k_tag
            save_plots_for_group(dfs, sub, k_tag, show=not args.no_show)

    save_plots_by_dt(frames, outdir, show=not args.no_show)



if __name__ == "__main__":
    main()
