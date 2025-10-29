# Morse Wavefunction Propagator

Symulator C++17 do jednowymiarowej ewolucji funkcji falowej w potencjale Morse’a. Program buduje trójdiagonalny hamiltonian na równomiernej siatce, diagonalizuje go spektralnie i propaguje znormalizowaną paczkę Gaussa metodą Taylora, RK4 lub Czebyszewa, zapisując diagnostykę do plików CSV.

## Spis treści

- [Przegląd](#przegląd)
- [Struktura projektu](#struktura-projektu)
- [Budowanie](#budowanie)
  - [Linux](#linux)
  - [Windows](#windows)
- [Uruchamianie symulacji](#uruchamianie-symulacji)
  - [Przykładowe polecenia](#przykładowe-polecenia)
  - [Skrypty PowerShell](#skrypty-powershell)
- [Parametry konfiguracyjne](#parametry-konfiguracyjne)
- [Wyjścia i logowanie](#wyjścia-i-logowanie)
- [Licencja](#licencja)

## Przegląd

Plik wykonywalny `morse` realizuje następujący przepływ pracy:

1. Odczytuje flagi CLI do struktury `Params` i tworzy symetryczną siatkę `Grid(N, xmax)` o kroku `dx = 2·xmax/(N-1)`.
2. Składa potencjał Morse’a `U(x)` oraz hamiltonian różnic skończonych `H = -½∂²/∂x² + U(x)` jako rzeczywistą macierz symetryczną, po czym przeskalowuje wektory własne do jednostkowej normy `L²` na siatce.
3. Diagonalizuje `H`, uzyskując pary własne, współczynniki spektralne i energię stanu podstawowego – dane te zasilają diagnostykę metryk θ w trakcie symulacji.
4. Inicjalizuje znormalizowaną rzeczywistą lub zespoloną paczkę Gaussa na wewnętrznych punktach siatki oraz skaluje hamiltonian ewolucji tak, aby spełnić `max(U)=Umax`, jeśli użytkownik tego żąda.
5. Wybiera ewolver (`TaylorEvolver`, `Rk4Evolver` lub `ChebyshevEvolver`), propaguje stan do czasu `tmax` z krokiem `dt`, a następnie zapisuje próbki CSV, opcjonalne zrzuty „wide” oraz podsumowania na konsoli.

Statyczne dane Morse’a (potencjał, wartości własne, wektory własne) są zapisywane przed ewolucją, chyba że ustawiono `--evolve_only`.
## Struktura projektu

```
.
├── CMakeLists.txt          # konfiguracja CMake (Eigen, flagi Release)
├── include/                # nagłówki publiczne
│   ├── core/               # narzędzia matematyczne, metryki spektralne, operacje trójdiagonalne
│   ├── evolve/             # interfejsy ewolverów i deklaracje integratorów
│   ├── cli.hpp             # struktura Params i parser CLI
│   ├── grid.hpp            # opis siatki jednorodnej
│   ├── hamiltonian.hpp     # pomocnicy do budowy hamiltonianu
│   ├── initial.hpp         # warunki początkowe Gaussa
│   └── ...
├── src/                    # implementacje odpowiadające nagłówkom
│   ├── core/               # normy, metryki θ, pomocnicze CSV
│   ├── evolve/             # integratory Taylora, RK4 i Czebyszewa oraz fabryka
│   ├── runtime_evolution.cpp  # skalowanie, logowanie, zrzuty wide
│   ├── main.cpp            # punkt wejścia CLI
│   └── ...
├── scripts/                # wykresy i wsadowe skrypty PowerShell
└── third_party/eigen-5.0.0 # opcjonalna kopia nagłówków Eigen
```

Najważniejsze moduły:

- `core/math_utils.*` – iloczyny skalarne, normy L², przekroje prawdopodobieństwa, mnożenie trójdiagonalne oraz ostrzeżenia o stabilności kroku czasowego.
- `core/io_utils.*` – obliczanie metryk θ, formatowanie wierszy CSV, logowanie konsolowe oraz generatory szerokich plików CSV.
- `core/spectral.*` – diagonalizacja macierzy i ewidencja współczynników spektralnych dla diagnostyki.
- `evolve/` – konkretne ewolwery: rozwinięcie Taylora z buforami roboczymi, Runge–Kutta 4. rzędu i adaptacyjny Czebyszew raportujący liczbę operacji macierz-wektor oraz opcjonalne `K_used`/`bn_ratio`.
- `hamiltonian.*` – budowa trójdiagonalnego operatora kinetycznego + potencjału oraz renormalizacja wektorów własnych do wag kwadratury siatki.
- `initial.*` – znormalizowane realne i zespolone paczki Gaussa na wewnętrznych węzłach siatki z opcjonalnym czynnikiem falowymi płaszczyzny.
- `io.*` – zapisy CSV, tabele widmowe i obsługa zrzutów wide (`*_abs2_wide.csv`, `*_re_wide.csv`, `*_im_wide.csv`).
- `runtime_evolution.*` – skalowanie potencjału, konstrukcja hamiltonianu, podłączanie ewolverów i zapis logów na dysk.

## Budowanie

Projekt korzysta z CMake, który wyszukuje bibliotekę Eigen (z pakietów systemowych lub dołączonego katalogu `third_party/eigen-5.0.0`). Flagi Release (`-O3`, `-DNDEBUG`, `EIGEN_NO_DEBUG`) są aktywne domyślnie.

### Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Plik wykonywalny znajdziesz w `build/morse`. Opcjonalnie uruchom test kontrolny: `ctest --test-dir build`.

### Windows

Skorzystaj z CMake wraz z MSVC lub MinGW (przykład w PowerShell):

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_CXX_FLAGS_RELEASE="/O2 /arch:AVX2 /DNDEBUG /D EIGEN_NO_DEBUG"
cmake --build . --config Release -j8
```

Wynikowy program to `build\Release\morse.exe` (lub `build\morse.exe` dla generatorów jednokonfiguracyjnych). Projekt ustawia `/W4 /O2` dla MSVC i korzysta wyłącznie z nagłówków Eigen, więc nie potrzeba dodatkowych bibliotek.

## Uruchamianie symulacji

Domyślnie wyniki trafiają do `results/` obok katalogu projektu. Ścieżka wyjściowa jest interpretowana względem pliku wykonywalnego, więc możesz podać zarówno ścieżki względne, jak i absolutne za pomocą `--outdir`.

### Przykładowe polecenia

Obliczanie samego widma (bez ewolucji w czasie):

```bash
./build/morse --N 2001 --xmax 30 --gamma 10 --first 10
```

Ewolucja czasowa metodą Taylora z zadanym rzędem:

```bash
./build/morse --evolve taylor --K 6 --dt 1e-5 --tmax 10
```

Ewolucja Czebyszewa z ostrą tolerancją:

```bash
./build/morse --evolve cheb --tol 1e-12 --dt 1e-4
```

Dodatkowe wskazówki:

- `--evolve_only` pomija ponowne generowanie statycznych plików Morse’a przy seryjnych uruchomieniach.
- `--quiet` wyłącza wypisywanie na konsolę przy zachowaniu plików CSV.
- `--csv results/custom.csv` nadpisuje domyślną nazwę pliku wynikowego.

### Skrypty PowerShell

Katalog `scripts/` zawiera wsadowe uruchamiacze dla przeglądów metod Taylora, RK4 i Czebyszewa. Możesz je uruchamiać w PowerShell 5.1+ na Windowsie oraz w PowerShell 7 (`pwsh`) na Linuksie/macOS; skrypty automatycznie wykrywają `morse.exe`, tworzą katalogi wyników dla każdej metody i startują zadania równoległe przy przypięciu zmiennych BLAS/OpenMP do jednego wątku.
Przykład (PowerShell w Windows):

```powershell 7
pwsh -File scripts/run_taylor.ps1
```

Na Linuksie zainstaluj PowerShell 7 (`sudo apt install powershell`) i uruchom identyczne polecenie `pwsh`. Przed dużymi seriami sprawdź każdą z wartości domyślnych (gamma, krok czasowy, logowanie) zapisanych w skryptach.


## Analiza i wykresy (Python)

Skrypt **`scripts/plot_compare_methods.py`** porównuje przebiegi metod **Taylor**, **RK4** i **Chebyshev** na podstawie logów CSV i generuje wykresy metryk **`norm_err`** oraz **`theta_abs`** (skala logarytmiczna).  
Dla Taylora automatycznie wybiera „najlepsze K” dla pary `(gamma, dt)` na podstawie minimalnego końcowego `theta_abs`.


### Wymagania

- Python **3.9+**
- Biblioteki: `pandas`, `matplotlib`, `numpy`

#### Instalacja (opcjonalnie wirtualne środowisko)

##### Linux / macOS
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -U pip pandas matplotlib numpy
```
##### Windows (PowerShell)
```bash
python -m venv .venv; .\.venv\Scripts\Activate.ps1
pip install -U pip pandas matplotlib numpy
```

#### Uruchomienie

Parametry:

- --root — katalog z logami (przeszukiwany rekurencyjnie; obsługa .csv/.tsv/*.log)
- --out — folder wyjściowy na wykresy (domyślnie ./plots_out)
- --dpi — rozdzielczość plików PNG (domyślnie 140)
- --verbose — dodatkowa diagnostyka


#### Przykłady

##### Linux / macOS
```bash
python scripts/plot_compare_methods.py --root results --out plots_out --dpi 160 --verbose
```

##### Windows (PowerShell):
```bash
python .\scripts\plot_compare_methods.py --root .\results --out .\plots_out --dpi 160 --verbose
```

Wykresy PNG zapisywane są w katalogu **plots_out/**

Nazwy plików:

- cmp_gamma<g>_dt<dt>_norm_err.png
- cmp_gamma<g>_dt<dt>_theta_abs.png


Konsola wypisuje krótkie podsumowania, np.:

>[INFO] gamma=10 dt=1e-4: Taylor K=6 (final theta_abs=8.75e-13, mean matvecs=4.00)


## Parametry konfiguracyjne

| Flaga | Domyślnie | Opis |
| --- | --- | --- |
| `--N` | `2001` | Liczba punktów siatki (z węzłami brzegowymi). |
| `--xmax` | `30.0` | Półszerokość obszaru `[-xmax, xmax]`.【|
| `--gamma` | `10.0` | Parametr potencjału Morse’a kontrolujący głębokość studni. |
| `--Umax` / `--Vcap` | `0.1` | Skalowanie potencjału używanego w ewolucji tak, aby `max(U)` nie przekraczało limitu; analiza widmowa pozostaje bez zmian.|
| `--first` | `10` | Liczba stanów własnych i energii zapisywanych do CSV oraz wypisywanych na stdout. |
| `--evolve <method>` | `taylor` | Aktywuje propagację w czasie i wybiera ewolver `taylor`, `rk4` lub `cheb` (niewrażliwe na wielkość liter). |
| `--dt` | `1e-5` | Stały krok czasowy dla wszystkich ewolverów. |
| `--tmax` | `10.0` | Łączny czas fizyczny; program wykonuje `round(tmax/dt)` kroków. |
| `--K` | `4` | Rząd obcięcia Taylora; przekazywany także do RK4 (ignorowany) oraz Czebyszewa (opcjonalny limit). |
| `--tol` | `1e-12` | Tolerancja metody Czebyszewa kontrolująca adaptacyjny stopień wielomianu. |
| `--log` / `--log-every` | `10000` | Zapis diagnostyki do konsoli i CSV co `log_every` kroków (zawsze pierwszy i ostatni). |
| `--csv-every` | `1` | Pomijanie wierszy CSV między logami – zapis co `csv_every` kroków.【F:include/cli.hpp†L17-L19】 |
| `--aggregate` | `false` | Agregacja czasu ściennego i liczby mnożeń w oknie logowania; w przeciwnym razie raportowany jest ostatni krok.|
| `--flush-every` | `1000` | Wymuszanie opróżnienia bufora CSV co podaną liczbę wierszy (0 wyłącza). |
| `--no-theta` | `false` | Pominięcie obliczeń metryk θ (tańsze, gdy diagnostyka spektralna nie jest potrzebna). |
| `--outdir` | `results` | Katalog wyjściowy; ścieżki względne są odnoszone do katalogu projektu. |
| `--csv <name>` | *(puste)* | Ręczne ustawienie nazwy pliku CSV; w przeciwnym razie generowany jest stem zależny od metody, `K`, `dt`, `gamma`, `N` i `xmax`.|
| `--evolve_only` | `false` | Pomija eksport statycznych plików Morse’a (`morse_*`) i uruchamia wyłącznie ewolucję.|
| `--quiet` | `false` | Wyłącza logowanie na konsolę przy zachowaniu plików CSV i wide. |
| `--wide` | `false` | Aktywuje szerokie zrzuty `|ψ|²` na wewnętrznej siatce (`*_abs2_wide.csv`). |
| `--wide-re` | `false` | Zapisuje dodatkowo część rzeczywistą (`*_re_wide.csv`). |
| `--wide-im` | `false` | Zapisuje dodatkowo część urojoną (`*_im_wide.csv`). |

## Wyjścia i logowanie

- Pliki statyczne: `morse_potential.csv`, `morse_eigenstates.csv`, `morse_energies_{num,anal}.csv` zapisywane przed ewolucją (chyba że ustawiono `--evolve_only`).
- Plik ewolucji (`*.csv`): kolumny `method, step, t, dt, dt_ms, matvecs, norm_err, theta_abs, theta_rel` oraz opcjonalnie `K_used` i `bn_ratio` dla metody Czebyszewa.
- Zrzuty wide: aktywowane flagą `--wide`. Plik `_abs2_wide.csv` zawiera `t` oraz gęstość prawdopodobieństwa na wewnętrznej siatce; opcje `--wide-re` i `--wide-im` tworzą dodatkowo `_re_wide.csv` i `_im_wide.csv`.
- Log konsolowy: formatowane wiersze z czasem trwania kroku, liczbą mnożeń macierz-wektor, dryftem `‖ψ‖₂` i metrykami θ dla każdego znacznika logowania (chyba że `--quiet`).

## Licencja

Projekt jest udostępniany na warunkach [licencji MIT](LICENSE).
