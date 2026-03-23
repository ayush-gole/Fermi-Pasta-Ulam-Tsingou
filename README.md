# Fermi-Pasta-Ulam-Tsingou (FPUT) Simulation

Numerical simulation of the FPUT problem using **Velocity Verlet integration**, written in C++ with Python for visualization of mode dynamics and recurrence.

Tracks how energy initially in a single normal mode evolves across all modes over time — and whether the system thermalizes or comes back to where it started.

---

## Outputs

| Plot | Description |
|------|-------------|
| ![Energy Evolution](Energy_Evolution_of_First_Three_modes.png) | Energy exchange between first three normal modes |
| ![Normalized Mode Energy](Normalized_Mode_Energy_Distribution.png) | Heatmap of normalized energy across all 32 modes |
| ![3D Energy Distribution](3D_Energy_Distribution.png) | 3D view of energy distribution across modes and time |
| ![R(t) Distribution](R_t__Distribution.png) | Recurrence function across all modes |
| ![3D Recurrence](3D_Recurance_Distribution.png) | 3D surface of the recurrence distribution |
| ![Recurrence vs Time](Recurance_Vs_Time_measure.png) | Scalar R(t) showing periodic return to initial state |

---

## Features

- Velocity Verlet integration for accurate time evolution
- Normal mode decomposition using discrete sine transform basis
- Tracks energy per mode and total energy over time
- Computes **recurrence function R(t)** — both scalar and per-mode
- Outputs full time series to CSV for easy post-processing
- Output filename changes automatically based on amplitude `A`

---

## How to Run

> Both files should be in the **same directory**

**Step 1 — Compile and run the C++ simulation:**
```bash
g++ -O2 -o FPUT FPUT.cpp
./FPUT
```
This generates a CSV named after the amplitude, for example:
- `FPUT_0.100000.csv` for A = 0.1
- `FPUT_1.000000.csv` for A = 1.0

**Step 2 — Update the file path in `plotting.py`:**

The script has a local path hardcoded, change this line to match your CSV:
```python
df = pd.read_csv('FPUT_0.100000.csv')
```

**Step 3 — Run the plotting script:**
```bash
python plotting.py
```

---

## Parameters

These are near the top of `FPUT.cpp`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 32 | Number of oscillators in the chain |
| `alpha` | 0.25 | Nonlinear coupling strength |
| `dt` | 0.05 | Integration time step |
| `maxT` | 50000 | Total simulation time |
| `A` | 0.1 | Initial amplitude (also sets output filename) |

---

## Tech Stack

- **C++** — Verlet integrator, normal mode analysis, CSV output
- **Python** — Data loading and visualization with `numpy`, `pandas`, `matplotlib`

---

## Physics Background

The FPUT problem is one of the earliest numerical experiments in physics (Fermi, Pasta, Ulam, Tsingou — 1955). A chain of oscillators with weak nonlinear coupling was expected to thermalize — energy starting in mode 1 should gradually spread across all modes. Instead the system kept returning close to its initial state, which was a pretty unexpected result at the time.

This behavior, now called FPUT recurrence, ended up being important for understanding solitons and the boundary between integrable and chaotic systems. The recurrence is measured using:

**R(t) = Σₖ [Qₖ(t)Qₖ(0) + Pₖ(t)Pₖ(0)] / Σₖ [Qₖ(0)² + Pₖ(0)²]**

R(t) = 1 means the system has fully returned to its initial state.
