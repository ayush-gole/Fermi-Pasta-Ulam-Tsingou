# Fermi-Pasta-Ulam-Tsingou (FPUT) Problem Simulation

A numerical simulation of the **FPUT problem** — one of the most historically significant problems in computational physics. Implemented in C++ using **Velocity Verlet integration**, with Python-based visualization of normal mode dynamics and recurrence phenomena.

The simulation reveals how energy, initially seeded into a single normal mode, fails to thermalize and instead **recurs periodically** — the famous FPUT recurrence paradox.

---

## Outputs

| Plot | Description |
|------|-------------|
| ![Energy Evolution](Energy_Evolution_of_First_Three_modes.png) | Energy exchange between first three normal modes over time |
| ![Normalized Mode Energy](Normalized_Mode_Energy_Distribution.png) | Heatmap of normalized energy across all 32 modes |
| ![3D Energy Distribution](3D_Energy_Distribution.png) | 3D surface of energy distribution across modes and time |
| ![R(t) Distribution](R_t__Distribution.png) | Recurrence function R(t) across all modes |
| ![3D Recurrence](3D_Recurance_Distribution.png) | 3D surface of recurrence distribution |
| ![Recurrence vs Time](Recurance_Vs_Time_measure.png) | Scalar recurrence R(t) showing periodic return to initial state |

---

## Features

- **Velocity Verlet integration** for accurate, symplectic time evolution
- Normal mode decomposition using **discrete sine transform (DST)** basis
- Tracks energy per normal mode Ek and total energy conservation
- Computes **recurrence function R(t)** — both scalar and per-mode
- Outputs full time series to CSV for flexible post-processing
- Configurable amplitude `A` (affects initial energy; CSV named accordingly)

---

## How to Run

> Both files must be in the **same directory**.

**Step 1 — Compile and run the C++ simulation:**
```bash
g++ -O2 -o FPUT FPUT.cpp
./FPUT
```
This generates a CSV file named after the amplitude, e.g.:
- `FPUT_0.100000.csv` for A = 0.1
- `FPUT_1.000000.csv` for A = 1.0

**Step 2 — Update the file path in `plotting.py`:**

Change this line to point to your generated CSV:
```python
df = pd.read_csv('FPUT_1.000000.csv')   # or FPUT_0.100000.csv
```

**Step 3 — Run the plotting script:**
```bash
python plotting.py
```

---

## Parameters (in `FPUT.cpp`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `N` | 32 | Number of oscillators in the chain |
| `alpha` | 0.25 | Nonlinear coupling strength |
| `dt` | 0.05 | Integration time step |
| `maxT` | 50000 | Total simulation time |
| `A` | 0.1 | Initial mode amplitude (also sets output filename) |

---

## Tech Stack

- **C++** — Velocity Verlet integrator, normal mode analysis, CSV output
- **Python** — Data loading and visualization (`numpy`, `pandas`, `matplotlib`)

---

## Physics Background

The **FPUT problem** was one of the first numerical experiments in physics (Fermi, Pasta, Ulam, Tsingou, 1955). A chain of *N* oscillators with a weak **nonlinear (quadratic) coupling**:

**aᵢ = (qᵢ₊₁ − 2qᵢ + qᵢ₋₁) + α[(qᵢ₊₁ − qᵢ)² − (qᵢ − qᵢ₋₁)²]**

Energy is initialized entirely in the **first normal mode** (k=1). Classical statistical mechanics predicts energy should eventually **equipartition** across all modes (thermalize). Instead, the system exhibits **quasi-periodic recurrence** — energy flows to a few low modes and returns nearly completely to the initial state, never thermalizing.

This unexpected result laid the groundwork for the discovery of **solitons** and the study of **chaos and integrability** in nonlinear systems.

The **recurrence function R(t)** quantifies how close the system returns to its initial state in phase space:

**R(t) = Σₖ [Qₖ(t)Qₖ(0) + Pₖ(t)Pₖ(0)] / Σₖ [Qₖ(0)² + Pₖ(0)²]**

R(t) = 1 means perfect recurrence to the initial configuration.
