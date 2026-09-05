# CHMRS-Glacial-Collapse
A Slow-Fast Dynamical Scaffold and Navier-Stokes Solver for Structural Glacial Range Failures.

# Coupled Hydro-Mechanical Retrogressive Slip (CHMRS) Framework
### A Multi-Paper Analytical and Numerical Trilogy for Structural Glacial Range Collapse

[![Language: Python](https://shields.io)](https://python.org)
[![Framework: LaTeX](https://shields.io)](https://latex-project.org)
[![License: MIT](https://shields.io)](LICENSE)

---

## 📌 Project Overview

This repository hosts the theoretical architecture, case study validations, and numerical simulation engines defining the **Coupled Hydro-Mechanical Retrogressive Slip (CHMRS)** framework. 

Traditional operational hazard tracking framework models treat high-altitude glacial mass movements through linear, precipitation-percolation pathways. CHMRS re-conceptualizes a heavily fractured glacier as an unattenuated, non-linear fluid-structure continuum. By coupling multi-dimensional incompressible **Navier-Stokes fluid dynamics** with the **Terzaghi-Glen Effective Pressure equations**, this project maps how episodic weather shock pulses trigger threshold crossings ($N \le 0$), forward-axis slurry venting, and backward-axis negative hydraulic shockwaves to drive cascading, retrogressive structural collapse.

---

## 📑 Repository Structure & Document Hierarchy

The research is organized as a unified, three-part academic series supported by an operational solver:

### 📄 [Paper 0] General Dynamical Scaffold
* **Title:** *A Slow-Fast Framework for Decadal-to-Millennial Rock-Ice Slope Collapse: Unifying Hydro-Mechanical and Thermal-Mechanical Failure as Fold-Bifurcation Instances*
* **Core Contribution:** Establishes the overarching meta-theoretical foundation. It models catastrophic slope failures as a slow-fast dynamical system where an internal degradation state variable ($D \in [0,1]$) is slaved to a fast, quasi-static force balance. It expands **Critical Slowing Down (CSD)** diagnostics from short-term seismic testing up to decadal-to-millennial scales.

### 📄 [Paper 1] Theoretical Architecture
* **Title:** *Coupled Hydro-Mechanical Retrogressive Slip (CHMRS): A Falsifiable Framework for Water-Pressure-Mediated Catastrophic Glacier Collapse*
* **Core Contribution:** Formulates the specific hydro-mechanical equations governing subglacial cavity over-pressurization. Critically, it outlines **five independently checkable Activation Criteria (C1–C5)** and strict precursor signature matrices, ensuring the hypothesis is independently checkable *before* an event is analyzed to eliminate circular configuration fits.

### 📄 [Paper 2] Case Study & Falsification Test
* **Title:** *Langtang Lirung, August 2026: An Evidence-First Test of the CHMRS Activation Criteria*
* **Core Contribution:** Executes an objective field-data test of the CHMRS criteria against the 26 August 2026 Langtang Lirung disaster. Because the failure plane consisted of high-grade metamorphic gneiss rather than erodible till (C2 fails), and zero rainfall occurred near the window (C5 fails), the paper demonstrates a **successful falsification**. It establishes that the event was a mechanically distinct thermal-mechanical permafrost-cemented rock-mass failure, vindicating the framework's mathematical boundaries.

### 💻 [Paper 2 / Code] Numerical Simulation Engine (`chmrs_simulation.py`)
* **Core Contribution:** An operational Python verification model that solves the 1D Navier-Stokes boundary equations under an episodic Heaviside shock loading function (e.g., a 50mm/2hr cloudburst) to dynamically map the **Time-to-Failure (TTF)** curve as $N \to 0$.

---

## 🧮 Mathematical Architecture

The underlying computational code tracks the fluid phase via an explicit finite-difference solution to the 1D incompressible Navier-Stokes equation containing advection, viscous diffusion, and gravity components:

$$\rho_w \left( \frac{\partial \vec{u}}{\partial t} + (\vec{u} \cdot \nabla)\vec{u} \right) = -\nabla P_w + \mu \nabla^2 \vec{u} + \rho_w \vec{g}$$

This fluid solver updates the Terzaghi effective pressure boundary layer across spatial grid nodes:

$$N = P_i - P_w \quad \text{where} \quad P_i = \rho_i g H$$

Flotation occurs dynamically, reducing the basal friction coefficient to zero, the instant:

$$N \le 0 \iff P_w \ge P_i$$

---

## 🚀 Execution & Quick Start (Paper 3 Solver)

### Prerequisites
The numerical script requires a standard Python environment with `numpy` and `matplotlib` installed.

```bash
pip install numpy matplotlib
```

### Running the Sensitivity Analysis
Execute the Python script from your terminal to simulate the subglacial pressure wave response to a Heaviside rainfall shock pulse:

```bash
python chmrs_simulation.py
```

### Expected Output
The script prints real-time step telemetry and generates a high-resolution analysis plot (`chmrs_simulation_output.png`) tracking two primary indicators:
1. **Basal Effective Pressure Drop:** Pinpoints the exact minute the flotation boundary threshold is breached.
2. **Terminal Outflow Velocity Pulse:** Captures the high-velocity forward venting phase prior to structural collapse.

---

## 📈 Early-Warning Application Matrix

For operational monitoring groups, the framework provides a qualitative diagnostic map derived from observed proglacial "surging pulses" in remote video footage:

* **Green Zone (Stable Drift):** Rhythmic, predictable water pulses separated by several hours. Internal conduits are venting efficiently.
* **Yellow Zone (Ambient Risk):** Pulse frequency compresses; river becomes heavily choked with grey silt. Frontal sapping is actively accelerating.
* **Orange Zone (Immediate Threat):** Downstream river flows suddenly drop or stop entirely. Internal conduits have suffered a roof collapse, causing subglacial water to back up under immense pressure head.
* **Red Zone (Critical Failure):** A sudden, hyper-dense wall of mud/debris advances downstream. Flotation limit has triggered a total retrogressive collapse.

---

## 👥 Citation & Academic Feedback

If you are reviewing this trilogy for submission or utilizing the CHMRS solver in secondary geomechanical modeling, please reference this unified repository:

```bibtex
@text{buchanan2026chmrs,
  author = {Buchanan, Paul},
  title = {Coupled Hydro-Mechanical Retrogressive Slip (CHMRS) Trilogy: A Unified Framework for Catastrophic Glacier Collapse},
  year = {2026},
  publisher = {GitHub},
  journal = {GitHub Repository},
  howpublished \(= {\url{https://github.com}} \)}
```
