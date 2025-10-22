# scripts/check_abs2_wide.py
import os, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

RESULTS = "results"
FIGS    = "figures"
EDGE_FRAC = 0.05   # доля точек с каждого края для оценки "утечки"

os.makedirs(FIGS, exist_ok=True)

def load_abs2(path):
    df = pd.read_csv(path, dtype=str, low_memory=False)
    t = pd.to_numeric(df['t'], errors='coerce').to_numpy()
    x_cols = df.columns[1:]
    x = np.array([float(c) for c in x_cols])
    Y = df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce').to_numpy()  # [nt, nx]
    return t, x, Y

def analyze_one(path):
    stem = os.path.basename(path).replace("_abs2_wide.csv", "")
    t, x, Y = load_abs2(path)
    dxs = np.diff(x)
    assert np.allclose(dxs, dxs[0], rtol=0, atol=1e-12), "x неравномерная сетка?"
    dx = dxs[0]

    # 1) Норма
    norms = Y.sum(axis=1) * dx
    drift = norms - 1.0
    print(f"[{stem}] norm: min={norms.min():.12f} max={norms.max():.12f}  max|drift|={np.max(np.abs(drift)):.3e}")

    # 2) Утечки на краях
    k = max(1, int(EDGE_FRAC * Y.shape[1]))
    edge_mass = (Y[:, :k].sum(axis=1) + Y[:, -k:].sum(axis=1)) * dx
    print(f"[{stem}] edge mass (≈leakage): max={edge_mass.max():.3e}")

    # 3) Моменты: x̄ и σ
    xbar = (Y @ x) * dx
    var  = (Y @ (x**2)) * dx - xbar**2
    var  = np.maximum(var, 0.0)
    sigma = np.sqrt(var)

    # --- графики ---
    def save_plot(y, ylabel, fname, yscale=None):
        plt.figure()
        plt.plot(t, y)
        plt.xlabel("t"); plt.ylabel(ylabel); plt.title(f"{ylabel} vs t — {stem}")
        if yscale: plt.yscale(yscale)
        plt.tight_layout()
        out = os.path.join(FIGS, f"{stem}_{fname}")
        plt.savefig(out, dpi=160, bbox_inches="tight"); plt.close()
        print("[saved]", out)

    save_plot(norms, "norm", "norm_vs_t_from_abs2.png")
    save_plot(drift, "norm-1", "norm_drift_vs_t_from_abs2.png")
    save_plot(edge_mass, "edge_mass", "edge_mass_vs_t.png", yscale="log")
    save_plot(xbar, "x̄(t)", "xbar_vs_t.png")
    save_plot(sigma, "σ(t)", "sigma_vs_t.png")

def main():
    files = sorted(glob.glob(os.path.join(RESULTS, "*_abs2_wide.csv")))
    if not files:
        print("[warn] no *_abs2_wide.csv in", RESULTS); return
    for f in files:
        try:
            analyze_one(f)
        except Exception as e:
            print("[error]", f, e)

if __name__ == "__main__":
    main()
