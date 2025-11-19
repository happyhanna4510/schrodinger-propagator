import re
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# ==============================
# USTAWIENIA
# ==============================

# 1) Ścieżka do logu theta-debug
log_path = Path(r"results\U0_-1\g10\dt1e-5\theta_debug_rk4_dt1e-5_g10_U0_-1.log")

# 2) Katalog główny na wykresy
plots_root = Path("plots_out")

# 3) Czy pokazywać okna matplotlib
show_plots = True

# ==============================
# ROZBIJANIE ŚCIEŻKI: U0 / gamma / dt / metoda / tag
# ==============================

try:
    dt_name = log_path.parent.name               # np. dt1e-5
    gamma_name = log_path.parent.parent.name     # np. g10
    u0_name = log_path.parent.parent.parent.name # np. U0_-100
except IndexError:
    print("[ERROR] Nie udalo sie wyciagnac U0/gamma/dt z log_path")
    sys.exit(1)

# metoda z nazwy pliku: theta_debug_<method>_...
m_method = re.search(r"theta_debug_([a-zA-Z0-9]+)_", log_path.name)
if m_method:
    method_name = m_method.group(1)  # cheb / rk4 / taylor
else:
    method_name = "unknown"

# tag dla podkatalogu (nazwa pliku bez rozszerzenia)
run_tag = log_path.stem  # np. theta_debug_cheb_dt1e-5_g10_U0_-100

# Katalog na wykresy:
# plots_out / U0_-100 / g10 / dt1e-5 / theta_debug_cheb_dt1e-5_g10_U0_-100
out_dir = plots_root / u0_name / gamma_name / dt_name / run_tag
out_dir.mkdir(parents=True, exist_ok=True)
print(f"[INFO] Katalog wyjsciowy: {out_dir}")

# ==============================
# PARSOWANIE LOGU
# ==============================

time_entries = {
    "step": [],
    "t_num": [],
    "t_ref": [],
    "theta_abs": [],
}

overlap_entries = {
    "t": [],
    "Re": [],
    "Im": [],
    "z_abs": [],
    "theta_raw": [],
    "theta_abs": [],
}

time_re = re.compile(
    r"step=([0-9]+).*t_num=([0-9eE\+\-\.]+).*t_ref=([0-9eE\+\-\.]+).*theta_abs=([0-9eE\+\-\.]+)"
)
overlap_re = re.compile(
    r"t=([0-9eE\+\-\.]+)\s+Re=([0-9eE\+\-\.]+)\s+Im=([0-9eE\+\-\.]+)\s+\|z\|=([0-9eE\+\-\.]+)\s+theta_raw=([0-9eE\+\-\.]+)\s+theta_abs=([0-9eE\+\-\.]+)"
)

with open(log_path, "r", encoding="utf-16") as f:
    for line in f:
        line = line.strip()
        if line.startswith("# [theta-debug-time]"):
            m = time_re.search(line)
            if m:
                step, t_num, t_ref, theta_abs = m.groups()
                time_entries["step"].append(int(step))
                time_entries["t_num"].append(float(t_num))
                time_entries["t_ref"].append(float(t_ref))
                time_entries["theta_abs"].append(float(theta_abs))
        elif line.startswith("# [theta-debug-overlap]"):
            m = overlap_re.search(line)
            if m:
                t, Re_z, Im_z, z_abs, theta_raw, theta_abs = m.groups()
                overlap_entries["t"].append(float(t))
                overlap_entries["Re"].append(float(Re_z))
                overlap_entries["Im"].append(float(Im_z))
                overlap_entries["z_abs"].append(float(z_abs))
                overlap_entries["theta_raw"].append(float(theta_raw))
                overlap_entries["theta_abs"].append(float(theta_abs))

print(f"Odczytano {len(time_entries['step'])} linii theta-debug-time")
print(f"Odczytano {len(overlap_entries['t'])} linii theta-debug-overlap")

if len(time_entries["step"]) == 0:
    print("\n[ERROR] Brak linii theta-debug-time w logu.")
    print("  -> Sprawdz, czy przy uruchomieniu byl ustawiony THETA_DEBUG=1.")
    print("  -> Sprawdz, czy log_path wskazuje na poprawny plik .log.")
    sys.exit(1)

# ==============================
# NUMPY
# ==============================

t_num = np.array(time_entries["t_num"])
t_ref = np.array(time_entries["t_ref"])
theta_abs_time = np.array(time_entries["theta_abs"])

dt_diff = t_num - t_ref

print("\nZgodnosc czasu:")
print("  max |t_num - t_ref| =", np.max(np.abs(dt_diff)))
print("  srednia |t_num - t_ref| =", np.mean(np.abs(dt_diff)))

# ==============================
# WYKRES 1: zgodnosc czasu
# ==============================

fig1 = plt.figure()
plt.plot(t_num, dt_diff, ".-")
plt.xlabel("t")
plt.ylabel("t_num - t_ref")
plt.title(f"Sprawdzenie zgodnosci czasu ({method_name}, {u0_name}, {gamma_name}, {dt_name})")
plt.grid(True)

fig1_path = out_dir / "wyrownanie_czasu.png"
fig1.savefig(fig1_path, dpi=150, bbox_inches="tight")
print(f"[INFO] Zapisano: {fig1_path}")

# ==============================
# WYKRESY 2 i 3: theta i |z|
# ==============================

if len(overlap_entries["t"]) > 0:
    t_ov  = np.array(overlap_entries["t"])
    z_abs = np.array(overlap_entries["z_abs"])
    theta_raw = np.array(overlap_entries["theta_raw"])
    theta_abs_ov = np.array(overlap_entries["theta_abs"])

    # --- Theta ---
    fig2 = plt.figure()
    plt.plot(t_ov, theta_raw, ".-", label="theta_raw (surowa)")
    plt.plot(t_ov, theta_abs_ov, ".-", label="theta_abs (bezwzgledna)")
    plt.xlabel("t")
    plt.ylabel("theta")
    plt.title(f"Zachowanie fazy w okolicy minimum ({method_name}, {u0_name}, {gamma_name}, {dt_name})")
    plt.legend()
    plt.grid(True)

    fig2_path = out_dir / "theta_okolica_minimum.png"
    fig2.savefig(fig2_path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Zapisano: {fig2_path}")

    # --- |z| ---
    fig3 = plt.figure()
    plt.plot(t_ov, z_abs, ".-")
    plt.xlabel("t")
    plt.ylabel("|z|")
    plt.title(f"Modul iloczynu skalarnego |<psi_ref|psi_num>| ({method_name}, {u0_name}, {gamma_name}, {dt_name})")
    plt.grid(True)

    fig3_path = out_dir / "modul_iloczynu_skal.png"
    fig3.savefig(fig3_path, dpi=150, bbox_inches="tight")
    print(f"[INFO] Zapisano: {fig3_path}")
else:
    print("\n[WARNING] Brak linii theta-debug-overlap.")
    print("  -> Byc moze THETA_DEBUG_WINDOW nie obejmuje danego zakresu czasu.")
    print("  -> W takim przypadku zapisujemy tylko wykres zgodnosci czasu.")

# ==============================
# Pokazywanie okien
# ==============================

if show_plots:
    plt.show()
else:
    plt.close("all")
