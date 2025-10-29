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
- [Wskazówki numeryczne](#wskazówki-numeryczne)
- [Licencja](#licencja)

## Przegląd

Plik wykonywalny `morse` realizuje następujący przepływ pracy:

1. Odczytuje flagi CLI do struktury `Params` i tworzy symetryczną siatkę `Grid(N, xmax)` o kroku `dx = 2·xmax/(N-1)`.【F:include/cli.hpp†L4-L35】【F:src/grid.cpp†L1-L4】
2. Składa potencjał Morse’a `U(x)` oraz hamiltonian różnic skończonych `H = -½∂²/∂x² + U(x)` jako rzeczywistą macierz symetryczną, po czym przeskalowuje wektory własne do jednostkowej normy `L²` na siatce.【F:src/morse_potential.cpp†L1-L112】【F:src/hamiltonian.cpp†L1-L34】
3. Diagonalizuje `H`, uzyskując pary własne, współczynniki spektralne i energię stanu podstawowego – dane te zasilają diagnostykę metryk θ w trakcie symulacji.【F:src/core/spectral.cpp†L1-L35】
4. Inicjalizuje znormalizowaną rzeczywistą lub zespoloną paczkę Gaussa na wewnętrznych punktach siatki oraz skaluje hamiltonian ewolucji tak, aby spełnić `max(U)=Umax`, jeśli użytkownik tego żąda.【F:src/initial.cpp†L1-L33】【F:src/runtime_evolution.cpp†L23-L66】
5. Wybiera ewolver (`TaylorEvolver`, `Rk4Evolver` lub `ChebyshevEvolver`), propaguje stan do czasu `tmax` z krokiem `dt`, a następnie zapisuje próbki CSV, opcjonalne zrzuty „wide” oraz podsumowania na konsoli.【F:src/evolve/evolve_factory.cpp†L66-L200】

Statyczne dane Morse’a (potencjał, wartości własne, wektory własne) są zapisywane przed ewolucją, chyba że ustawiono `--evolve_only`.【F:src/main.cpp†L30-L52】【F:src/morse_static.cpp†L12-L71】

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

- `core/math_utils.*` – iloczyny skalarne, normy L², przekroje prawdopodobieństwa, mnożenie trójdiagonalne oraz ostrzeżenia o stabilności kroku czasowego.【F:src/core/math_utils.cpp†L1-L76】
- `core/io_utils.*` – obliczanie metryk θ, formatowanie wierszy CSV, logowanie konsolowe oraz generatory szerokich plików CSV.【F:src/core/io_utils.cpp†L1-L120】
- `core/spectral.*` – diagonalizacja macierzy i ewidencja współczynników spektralnych dla diagnostyki.【F:src/core/spectral.cpp†L1-L35】
- `evolve/` – konkretne ewolwery: rozwinięcie Taylora z buforami roboczymi, Runge–Kutta 4. rzędu i adaptacyjny Czebyszew raportujący liczbę operacji macierz-wektor oraz opcjonalne `K_used`/`bn_ratio`.【F:include/evolve/evolver_base.hpp†L5-L40】【F:src/evolve/evolve_factory.cpp†L118-L193】
- `hamiltonian.*` – budowa trójdiagonalnego operatora kinetycznego + potencjału oraz renormalizacja wektorów własnych do wag kwadratury siatki.【F:src/hamiltonian.cpp†L1-L34】
- `initial.*` – znormalizowane realne i zespolone paczki Gaussa na wewnętrznych węzłach siatki z opcjonalnym czynnikiem falowymi płaszczyzny.【F:src/initial.cpp†L1-L33】
- `io.*` – zapisy CSV, tabele widmowe i obsługa zrzutów wide (`*_abs2_wide.csv`, `*_re_wide.csv`, `*_im_wide.csv`).【F:src/io.cpp†L1-L205】
- `runtime_evolution.*` – skalowanie potencjału, konstrukcja hamiltonianu, podłączanie ewolverów i zapis logów na dysk.【F:src/runtime_evolution.cpp†L1-L108】

## Budowanie

Projekt korzysta z CMake, który wyszukuje bibliotekę Eigen (z pakietów systemowych lub dołączonego katalogu `third_party/eigen-5.0.0`). Flagi Release (`-O3`, `-DNDEBUG`, `EIGEN_NO_DEBUG`) są aktywne domyślnie.【F:CMakeLists.txt†L1-L63】

### Linux

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Plik wykonywalny znajdziesz w `build/morse`. Opcjonalnie uruchom test kontrolny: `ctest --test-dir build`.【F:CMakeLists.txt†L65-L78】

### Windows

Skorzystaj z CMake wraz z MSVC lub MinGW (przykład w PowerShell):

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Wynikowy program to `build\Release\morse.exe` (lub `build\morse.exe` dla generatorów jednokonfiguracyjnych). Projekt ustawia `/W4 /O2` dla MSVC i korzysta wyłącznie z nagłówków Eigen, więc nie potrzeba dodatkowych bibliotek.【F:CMakeLists.txt†L23-L55】

## Uruchamianie symulacji

Domyślnie wyniki trafiają do `results/` obok katalogu projektu. Ścieżka wyjściowa jest interpretowana względem pliku wykonywalnego, więc możesz podać zarówno ścieżki względne, jak i absolutne za pomocą `--outdir`.【F:src/main.cpp†L32-L46】【F:include/paths.hpp†L7-L15】

### Przykładowe polecenia

Obliczanie samego widma (bez ewolucji w czasie):

```bash
./build/morse --N 2001 --xmax 30 --gamma 10 --first 10
```

Ewolucja czasowa metodą Taylora z zadanym rzędem:

```bash
./build/morse --evolve taylor --K 6 --dt 1e-5 --tmax 1.0
```

Ewolucja Czebyszewa z ostrą tolerancją:

```bash
./build/morse --evolve cheb --tol 1e-12 --dt 1e-4 --tmax 1.0
```

Dodatkowe wskazówki:

- `--evolve_only` pomija ponowne generowanie statycznych plików Morse’a przy seryjnych uruchomieniach.【F:src/main.cpp†L42-L48】
- `--quiet` wyłącza wypisywanie na konsolę przy zachowaniu plików CSV.【F:include/cli.hpp†L28-L31】
- `--csv results/custom.csv` nadpisuje domyślną nazwę pliku wynikowego.【F:src/cli.cpp†L45-L50】

### Skrypty PowerShell

Katalog `scripts/` zawiera wsadowe uruchamiacze dla przeglądów metod Taylora, RK4 i Czebyszewa. Możesz je uruchamiać w PowerShell 5.1+ na Windowsie oraz w PowerShell 7 (`pwsh`) na Linuksie/macOS; skrypty automatycznie wykrywają `morse.exe`, tworzą katalogi wyników dla każdej metody i startują zadania równoległe przy przypięciu zmiennych BLAS/OpenMP do jednego wątku.【F:scripts/run_taylor.ps1†L1-L146】【F:scripts/run_cheb.ps1†L1-L49】

Przykład (PowerShell w Windows):

```powershell
pwsh -File scripts/run_taylor.ps1
```

Na Linuksie zainstaluj PowerShell (`sudo apt install powershell`) i uruchom identyczne polecenie `pwsh`. Przed dużymi seriami sprawdź każdą z wartości domyślnych (gamma, krok czasowy, logowanie) zapisanych w skryptach.【F:scripts/run_taylor.ps1†L24-L94】

## Parametry konfiguracyjne

| Flaga | Domyślnie | Opis |
| --- | --- | --- |
| `--N` | `2001` | Liczba punktów siatki (z węzłami brzegowymi).【F:include/cli.hpp†L5-L6】 |
| `--xmax` | `30.0` | Półszerokość obszaru `[-xmax, xmax]`.【F:include/cli.hpp†L6-L7】 |
| `--gamma` | `10.0` | Parametr potencjału Morse’a kontrolujący głębokość studni.【F:include/cli.hpp†L7-L8】 |
| `--Umax` / `--Vcap` | `0.1` | Skalowanie potencjału używanego w ewolucji tak, aby `max(U)` nie przekraczało limitu; analiza widmowa pozostaje bez zmian.【F:include/cli.hpp†L8-L9】【F:src/runtime_evolution.cpp†L26-L38】 |
| `--first` | `10` | Liczba stanów własnych i energii zapisywanych do CSV oraz wypisywanych na stdout.【F:include/cli.hpp†L21-L22】【F:src/morse_static.cpp†L32-L53】 |
| `--evolve <method>` | `taylor` | Aktywuje propagację w czasie i wybiera ewolver `taylor`, `rk4` lub `cheb` (niewrażliwe na wielkość liter).【F:src/cli.cpp†L34-L54】【F:src/evolve/evolve_factory.cpp†L102-L152】 |
| `--dt` | `1e-5` | Stały krok czasowy dla wszystkich ewolverów.【F:include/cli.hpp†L13-L14】 |
| `--tmax` | `10.0` | Łączny czas fizyczny; program wykonuje `round(tmax/dt)` kroków.【F:include/cli.hpp†L14-L15】【F:src/runtime_evolution.cpp†L57-L64】 |
| `--K` | `4` | Rząd obcięcia Taylora; przekazywany także do RK4 (ignorowany) oraz Czebyszewa (opcjonalny limit).【F:include/cli.hpp†L12-L13】【F:src/evolve/evolve_factory.cpp†L137-L155】 |
| `--tol` | `1e-12` | Tolerancja metody Czebyszewa kontrolująca adaptacyjny stopień wielomianu.【F:include/cli.hpp†L15-L16】【F:include/evolve/evolver_base.hpp†L9-L18】 |
| `--log` / `--log-every` | `10000` | Zapis diagnostyki do konsoli i CSV co `log_every` kroków (zawsze pierwszy i ostatni).【F:include/cli.hpp†L16-L20】【F:src/evolve/evolve_factory.cpp†L156-L200】 |
| `--csv-every` | `1` | Pomijanie wierszy CSV między logami – zapis co `csv_every` kroków.【F:include/cli.hpp†L17-L19】 |
| `--aggregate` | `false` | Agregacja czasu ściennego i liczby mnożeń w oknie logowania; w przeciwnym razie raportowany jest ostatni krok.【F:include/cli.hpp†L18-L19】【F:src/evolve/evolve_factory.cpp†L86-L115】 |
| `--flush-every` | `1000` | Wymuszanie opróżnienia bufora CSV co podaną liczbę wierszy (0 wyłącza).【F:include/cli.hpp†L19-L20】【F:src/evolve/evolve_factory.cpp†L173-L188】 |
| `--no-theta` | `false` | Pominięcie obliczeń metryk θ (tańsze, gdy diagnostyka spektralna nie jest potrzebna).【F:include/cli.hpp†L20-L21】【F:src/evolve/evolve_factory.cpp†L166-L191】 |
| `--outdir` | `results` | Katalog wyjściowy; ścieżki względne są odnoszone do katalogu projektu.【F:include/cli.hpp†L22-L24】【F:include/paths.hpp†L7-L15】 |
| `--csv <name>` | *(puste)* | Ręczne ustawienie nazwy pliku CSV; w przeciwnym razie generowany jest stem zależny od metody, `K`, `dt`, `gamma`, `N` i `xmax`.【F:include/cli.hpp†L23-L24】【F:src/io.cpp†L125-L168】 |
| `--evolve_only` | `false` | Pomija eksport statycznych plików Morse’a (`morse_*`) i uruchamia wyłącznie ewolucję.【F:include/cli.hpp†L26-L27】【F:src/main.cpp†L40-L48】 |
| `--quiet` | `false` | Wyłącza logowanie na konsolę przy zachowaniu plików CSV i wide.【F:include/cli.hpp†L27-L28】【F:src/evolve/evolve_factory.cpp†L184-L193】 |
| `--wide` | `false` | Aktywuje szerokie zrzuty `|ψ|²` na wewnętrznej siatce (`*_abs2_wide.csv`).【F:include/cli.hpp†L28-L30】【F:src/runtime_evolution.cpp†L68-L90】 |
| `--wide-re` | `false` | Zapisuje dodatkowo część rzeczywistą (`*_re_wide.csv`).【F:include/cli.hpp†L29-L31】【F:src/io.cpp†L169-L205】 |
| `--wide-im` | `false` | Zapisuje dodatkowo część urojoną (`*_im_wide.csv`).【F:include/cli.hpp†L29-L31】【F:src/io.cpp†L169-L205】 |

## Wyjścia i logowanie

- Pliki statyczne: `morse_potential.csv`, `morse_eigenstates.csv`, `morse_energies_{num,anal}.csv` zapisywane przed ewolucją (chyba że ustawiono `--evolve_only`).【F:src/morse_static.cpp†L44-L70】
- Plik ewolucji (`*.csv`): kolumny `method, step, t, dt, dt_ms, matvecs, norm_err, theta_abs, theta_rel` oraz opcjonalnie `K_used` i `bn_ratio` dla metody Czebyszewa.【F:src/core/io_utils.cpp†L91-L120】
- Zrzuty wide: aktywowane flagą `--wide`. Plik `_abs2_wide.csv` zawiera `t` oraz gęstość prawdopodobieństwa na wewnętrznej siatce; opcje `--wide-re` i `--wide-im` tworzą dodatkowo `_re_wide.csv` i `_im_wide.csv`.【F:src/io.cpp†L169-L205】
- Log konsolowy: formatowane wiersze z czasem trwania kroku, liczbą mnożeń macierz-wektor, dryftem `‖ψ‖₂` i metrykami θ dla każdego znacznika logowania (chyba że `--quiet`).【F:src/core/io_utils.cpp†L121-L155】

## Wskazówki numeryczne

- **Stabilność kroku czasowego** – Program ostrzega, gdy `dt · ΔE > 0.5`; w takiej sytuacji zmniejsz `dt` lub zwiększ rząd Taylora `K`, aby utrzymać stabilność rozwinięć Taylora i Czebyszewa.【F:src/core/math_utils.cpp†L67-L75】
- **Metryki θ** – `theta_abs` to kwadrat błędu bezwzględnego między współczynnikami spektralnymi a stanem po ewolucji numerycznej; `theta_rel` to wariant względny. Wyłącz `--no-theta`, aby pominąć projekcje na wektory własne.【F:src/core/io_utils.cpp†L32-L90】
- **Budowanie w Release** – Do długich symulacji kompiluj z optymalizacjami (`-O3`, `-DNDEBUG`); buildy debug znacznie spowalniają gęsty algorytm Eigena.【F:CMakeLists.txt†L1-L47】
- **Wielowątkowość** – Eigen domyślnie działa jednosesyjnie. Skrypty PowerShell przypinają zmienne BLAS/OpenMP do wartości 1, aby uniknąć nadsubskrypcji podczas równoległych uruchomień.【F:scripts/run_taylor.ps1†L56-L65】

## Licencja

Projekt jest udostępniany na warunkach [licencji MIT](LICENSE).【F:LICENSE†L1-L19】
