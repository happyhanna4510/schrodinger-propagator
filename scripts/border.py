import numpy as np
import pandas as pd

PATH = r"results\taylor_K6_dt1e-4_g10_N2001_x30_abs2_wide.csv"

df = pd.read_csv(PATH)
if df.empty:
    raise SystemExit("file is empty")

last = df.iloc[-1]
value_cols = df.columns[1:]
values = last[value_cols].astype(float).to_numpy()

k = 5
left = values[:k]
right = values[-k:]
print("left edges:", left)
print("right edges:", right)
print("max edge abs2:", float(np.max([left.max(), right.max()])))
