# run_energy_plots.ps1
# Рисует E(t) для:
# U0 = -30, -1, 1
# методы: rk4, taylor_K4, cheb

# --- U0 = -30 ---

python .\scripts\plot_energy.py `
  .\results\U0_-30\g10\dt1e-5\rk4_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_-30\g10\dt1e-5




python .\scripts\plot_energy.py `
  .\results\U0_-30\g10\dt1e-5\taylor_K4_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_-30\g10\dt1e-5

python .\scripts\plot_energy.py `
  .\results\U0_-30\g10\dt1e-5\cheb_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_-30\g10\dt1e-5


# --- U0 = -1 ---

python .\scripts\plot_energy.py `
  .\results\U0_-1\g10\dt1e-5\rk4_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_-1\g10\dt1e-5

python .\scripts\plot_energy.py `
  .\results\U0_-1\g10\dt1e-5\taylor_K4_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_-1\g10\dt1e-5

python .\scripts\plot_energy.py `
  .\results\U0_-1\g10\dt1e-5\cheb_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_-1\g10\dt1e-5


# --- U0 = 1 ---

python .\scripts\plot_energy.py `
  .\results\U0_1\g10\dt1e-5\rk4_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_1\g10\dt1e-5

python .\scripts\plot_energy.py `
  .\results\U0_1\g10\dt1e-5\taylor_K4_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_1\g10\dt1e-5

python .\scripts\plot_energy.py `
  .\results\U0_1\g10\dt1e-5\cheb_dt1e-5_g10_N2001_x20_energy.csv `
  --out .\plots_out\U0_1\g10\dt1e-5
