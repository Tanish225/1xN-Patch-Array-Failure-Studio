# ArrayFailureStudioN: 1 × N Patch Array Failure Analysis & Reconfiguration

A MATLAB GUI application designed for failure analysis, beam-preserving reconfiguration, and S-parameter evaluation of $1 \times N$ microstrip patch antenna arrays.

---

## Technical Overview

`ArrayFailureStudioN` models element failures in a linear antenna array and re-synthesizes active element excitation weights to compensate for gain loss and prevent main-beam steering drift. The algorithm formulates the active set optimization as a constrained weighted least-squares problem, solved in closed form using a Karush-Kuhn-Tucker (KKT) linear system.

### Mathematical Formulation

Given a total set of elements $N$, a set of faulted elements $F$, and a set of active healthy elements $S$ ($\vert{}S\vert{} \ge 2$), the excitation vector for active elements $w_S$ is solved via:

$$\min_{w_S} \left\Vert{} W^{1/2} (G_S w_S - b) \right\Vert{}^2$$

$$\text{subject to} \quad \begin{cases} \text{Keep peak value of ideal beam (hard constraint)} \\ \left. \frac{\partial F}{\partial \theta} \right\vert{}_{\text{peak}} = \left. \frac{\partial F}{\partial \phi} \right\vert{}_{\text{peak}} = 0 \quad (\text{soft/locked constraint}) \end{cases}$$

Where:

* $G_S$: Array manifold matrix corresponding to healthy elements $S$.
* $b = F_{\text{ideal}} - G_F w_F$: Residual target pattern subtracting contributions from faulted weights $w_F$.
* $W$: Solid-angle weighting diagonal matrix incorporating main-beam emphasis.

---

## System Requirements

* **MATLAB Version:** R2020a or newer (requires `uifigure` and App Designer UI component support).
* **Toolbox Requirements:** None (built entirely using standard MATLAB language features and native graphics).

---

## Quick Start

1. Place `ArrayFailureStudioN.m` in your MATLAB working directory or path.
2. Launch the interactive dashboard by executing the function in the Command Window:
```matlab
ArrayFailureStudioN

```



---

## Key Features

### 1. Active Set Search Optimization

* **Auto Mode:** Performs an exhaustive search across all valid healthy subsets when $\vert{}S\vert{} \le 8$; defaults to greedy backward elimination for larger active element counts.
* **Exhaustive Search:** Evaluates every subset combination to locate the global optimal fidelity configuration (hard limit up to 12 elements).
* **Greedy Backward Elimination:** Iteratively removes active elements whose exclusion minimizes pattern correlation degradation.

### 2. External Data Integration (CST Studio Suite)

* **Element Radiation Patterns:** Import custom 1D far-field cut files (`.txt`, `.csv`, `.dat`) for E-plane ($\phi = 0^\circ$) and H-plane ($\phi = 90^\circ$) cuts exported from CST Studio Suite, or fall back to an analytical $\cos^q(\theta)$ model.
* **Mutual Coupling & S-Parameters:** Load Touchstone files (`.sNp`) directly from CST to incorporate full mutual coupling matrices ($S_{ij}$) into reflection and active port power metrics. Alternatively, use the built-in single-port parametric resonator model.

### 3. Visualization & Diagnostics Suite

* **Live 3D & Geometry View:** Real-time 3D beam surface plots alongside physical substrate geometry visualization showing patch status (Healthy, Faulted, Switched Off).
* **Polar Comparison:** Direct overlay of ideal, failed, and compensated patterns in both the main beam plane ($\phi = \phi_0$) and the orthogonal plane ($\phi = \phi_0 + 90^\circ$).
* **Angle Inspector:** Interactive slider-driven angular sampling yielding local field values, gains ($\text{dBi}$), and phase profiles.
* **Beam Proof Module:** Computes beam peak theta values across all azimuth cuts ($\phi \in [0^\circ, 345^\circ]$) to verify directional stability under a $0.5^\circ$ tolerance threshold.

---

## Input Control Settings

| Category | Parameter | Description | Default | Range |
| --- | --- | --- | --- | --- |
| **Array & Substrate** | `Number of patches N` | Total array element count | `20` | $2 - 32$ |
|  | `Frequency f` | Operating center frequency ($\text{GHz}$) | `10` | $1 - 100\text{ GHz}$ |
|  | `Substrate eps_r` | Relative permittivity $\varepsilon_r$ | `4.4` | $1 - 20$ |
|  | `Substrate height h` | Substrate thickness ($\text{mm}$) | `1.6` | $0.1 - 10\text{ mm}$ |
|  | `Spacing d` | Center-to-center pitch ($\lambda_0$) | `0.5` | $0.2 - 2.0\ \lambda_0$ |
|  | `Steer theta0 / phi0` | Main-beam target steering angles ($\text{deg}$) | `30, 0` | $\theta \in [0, 90]^\circ, \phi \in [-360, 360]^\circ$ |
| **Element Pattern** | `Element source` | Selection between $\cos^q$ model or imported CST cuts | `cos^q model` | Enumerated |
|  | `cos^q exponent q` | Pattern roll-off exponent | `1.5` | $0 - 10$ |
| **Failure & Algorithm** | `Failed patches` | Selection list of defective array elements | None | Multiselect listbox |
|  | `Failure severity` | Defect fraction ($0 = \text{healthy}, 1 = \text{dead}$) | `1.0` | $0.0 - 1.0$ |
|  | `Beam-lock weight` | KKT constraint weight for steering preservation | `1000` | $0 - 10^{12}$ |
|  | `Reconfiguration search` | Active set selection heuristic | `Auto` | `Auto`, `Exhaustive`, `Greedy` |

---

## Data Export Capabilities

Clicking **Export results** generates a structured output package consisting of:

* `.mat`: Full workspace archive containing configuration settings, array matrices, and current metric states.
* `_pattern_3D.csv`: Tabulated 3D radiation dataset ($\theta, \phi$, Ideal dB, Failed dB, Compensated dB).
* `_failure_analysis.csv`: Batch output log of evaluated failure cases and fidelity ratings.
* `_beam_proof.csv`: Deviation log across all evaluated azimuth cuts.
* `_S11_current.csv`: Port reflection spectra across frequency.
* `_dashboard.png`: High-resolution render ($200\text{ DPI}$) of the application interface.
* `_report.txt`: Summary text report detailing array metrics, parameters, and complex excitation coefficients.
