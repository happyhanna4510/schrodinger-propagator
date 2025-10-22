# scripts/compare_abs2_pairs.py
import glob, os
import pandas as pd
import numpy as np

def load_abs2(p):
    df = pd.read_csv(p)
    t  = df['t'].to_numpy(float)
    Y  = df.drop(columns=['t']).to_numpy(float)  # [nt, nx]
    return t, Y

files = sorted(glob.glob("results/*_abs2_wide.csv"))
for i in range(len(files)):
    for j in range(i+1, len(files)):
        t1, Y1 = load_abs2(files[i])
        t2, Y2 = load_abs2(files[j])
        if len(t1)!=len(t2) or not np.allclose(t1, t2):
            print(f"[skip] different t-grid: {os.path.basename(files[i])} vs {os.path.basename(files[j])}")
            continue
        if Y1.shape != Y2.shape:
            print(f"[skip] different shapes")
            continue
        diff = np.max(np.abs(Y1 - Y2))
        print(f"{os.path.basename(files[i])} vs {os.path.basename(files[j])}: max |Δ| = {diff:.3e}")
