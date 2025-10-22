# scripts/plot_errors_by_stem.py
import glob, os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def load_log(p):
    df = pd.read_csv(p)
    # ожидаемые колонки: step,t,norm,err_true,theta (+возможно p0,err_phi0)
    return (df['t'].to_numpy(float),
            df['norm'].to_numpy(float),
            df['err_true'].to_numpy(float),
            df['theta'].to_numpy(float))

# собираем логи (без wide-файлов)
logs = sorted(glob.glob("results/*.csv"))
logs = [p for p in logs if not any(s in p for s in ["_abs2_wide","_re_wide","_im_wide"])]

series = []
for p in logs:
    try:
        t, normv, err_true, theta = load_log(p)
        stem = os.path.basename(p).replace(".csv","")
        series.append((stem, t, normv, err_true, theta))
    except Exception:
        pass

def plot_many(title, ylabel, extractor, outname):
    plt.figure()
    for stem, t, normv, err_true, theta in series:
        y = extractor(normv, err_true, theta)
        plt.plot(t, y, label=stem)
    plt.xlabel("t")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(fontsize=8)
    os.makedirs("figures", exist_ok=True)
    plt.savefig(os.path.join("figures", outname), dpi=160, bbox_inches="tight")
    plt.close()
    print("[saved]", outname)

# 0) Норма как есть
plot_many("Norm vs t", "norm", lambda n,e,th: n, "norm_vs_t.png")

# 1) Дрейф нормы (norm - 1)
plot_many("Norm drift vs t", "norm - 1", lambda n,e,th: n - 1.0, "norm_drift_vs_t.png")

# 2) err_true
plot_many("err_true vs t", "err_true", lambda n,e,th: e, "err_true_vs_t.png")

# 3) Theta
plot_many("Theta vs t", "Theta", lambda n,e,th: th, "theta_vs_t.png")
