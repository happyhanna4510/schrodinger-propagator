import pandas as pd
df = pd.read_csv(r"results\taylor_K6_dt1e-4_g10_N2001_x30_abs2_wide.csv")
last = df.iloc[-1]
# первые 5 и последние 5 узлов:
left = last[['y0','y1','y2','y3','y4']].values
right = last[[f'y{i}' for i in range(len(last)-1 - 4, len(last)-1)]].values  # -1 для 't'
print("left edges:", left)
print("right edges:", right)
print("max edge abs2:", max(left.max(), right.max()))
