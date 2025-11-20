# run_energies.ps1
# Последовательно считает энергии для:
# методы: rk4, taylor K=4, cheb
# U0: -30, -1, 1
# dt = 1e-5, gamma = 10, N = 2001, xmax = 20, tmax = 10
# Вывод полностью подавлен (> $null 2>&1)

# --- RK4 ---

& .\build\Release\morse.exe `
  --evolve rk4 `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 -30 `
  --outdir .\results\U0_-30\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1


& .\build\Release\morse.exe `
  --evolve rk4 `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 -1 `
  --outdir .\results\U0_-1\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1

& .\build\Release\morse.exe `
  --evolve rk4 `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 1 `
  --outdir .\results\U0_1\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1


# --- TAYLOR K=4 ---

& .\build\Release\morse.exe `
  --evolve taylor --K 4 `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 -30 `
  --outdir .\results\U0_-30\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1

& .\build\Release\morse.exe `
  --evolve taylor --K 4 `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 -1 `
  --outdir .\results\U0_-1\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1

& .\build\Release\morse.exe `
  --evolve taylor --K 4 `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 1 `
  --outdir .\results\U0_1\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1


# --- CHEB (если у тебя другой флаг, напр. --evolve chebyshev, замени здесь) ---

& .\build\Release\morse.exe `
  --evolve cheb `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 -30 `
  --outdir .\results\U0_-30\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1

& .\build\Release\morse.exe `
  --evolve cheb `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 -1 `
  --outdir .\results\U0_-1\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1

& .\build\Release\morse.exe `
  --evolve cheb `
  --gamma 10 --N 2001 --xmax 20 `
  --tmax 10 --dt 1e-5 `
  --init complex-gauss --k0 0 --U0 1 `
  --outdir .\results\U0_1\g10\dt1e-5 `
  --csv-every 0 `
  --log-energy `
  > $null 2>&1
