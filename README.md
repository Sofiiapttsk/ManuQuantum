# ManuQuantum (Quantum Crystal Growth Simulation (Earth vs Microgravity with Acoustic Levitation))

An educational, interactive Python / PyQt5 + OpenGL application that visualizes conceptual crystal growth for advanced quantum materials (e.g. Bi₂Se₃, YBCO, α‑RuCl₃) under:
- Earth gravity (buoyant convection + dopant advection)
- Microgravity (diffusion‑dominated transport)
- Acoustic levitation (containerless processing & stabilization)
> NOTE: This is a pedagogical visualization prototype, not a validated scientific CFD / phase‑field solver.

-

## Key Features

| Category | Feature |
|----------|---------|
| Environment Modes | Earth (1g) vs Space (µg) |
| Transport | 2D axisymmetric dopant diffusion + semi‑Lagrangian advection in Earth mode |
| Convection | Synthetic incompressible roll field (Earth) vs suppressed diffusion (Space) |
| Acoustic Field | Standing-wave inspired potential: surface modulation + convection damping |
| Dopant Incorporation | Interfacial depletion band with kinetic uptake rate |
| Defects | Stochastic generation (environment dependent spatial bias) |
| Crystal Quality Metric | Composite: dopant uniformity, temperature proximity, acoustic stabilization, defect probability |
| Visualization | 3D melt & crystal, flow lines, dopant markers, defects, acoustic sources |
| 2D Panel | Stylized convection/diffusion particle tracers and dopant proxies |
| Materials | Bi₂Se₃ (Topological Insulator), YBCO (High‑Tc), α‑RuCl₃ (Spin Liquid candidate) |

-

## Physical Model (Simplified)

| Aspect | Implementation | Not Modeled |
|--------|----------------|-------------|
| Fluid Flow | Analytic stream function → toroidal roll (Earth) | Navier–Stokes, buoyancy feedback, turbulence |
| Diffusion | Explicit Laplacian (axisymmetric r–z grid) | Temperature-dependent D, anisotropy |
| Advection | Semi‑Lagrangian backtrace (Earth mode) | Full mass conservation in 3D |
| Acoustic Field | Superposed damped spherical wave heuristic + Y₂₀/Y₂₂ surface perturbation | Helmholtz / full acoustic pressure solution |
| Interface | Radial growth fraction scales geometry; depletion band removes dopants | Faceting, curvature effects, Gibbs–Thomson |
| Defects | Stochastic spawn scaled by environment probability | Energetics, dislocation motion, clustering kinetics |
| Quality | Weighted scalar (uniformity + temp + acoustic + defect) | Mechanical stress, thermal shock, phase impurities |

-

## 📦 Requirements

| Dependency | Purpose |
|------------|---------|
| Python 3.9+ | Core runtime |
| PyQt5 | GUI framework |
| PyOpenGL | Rendering API |
| numpy | Numeric grids & statistics |
| OpenGL driver | Hardware acceleration (legacy immediate mode) |

### Install

```bash
python -m venv venv
source venv/bin/activate          # (Linux/macOS)
# venv\Scripts\activate           # (Windows PowerShell)

pip install --upgrade pip
pip install PyQt5 PyOpenGL numpy
```

### Run

```bash
python quantumcrystalgrowthearthspaceEnhanced.py
```

If Unicode issues appear (e.g. misrendered α, ₂):
```bash
export PYTHONUTF8=1
```

-

## 🖥️ User Interface Guide

### Tabs

| Tab | Purpose |
|-----|---------|
| 3D Crystal Growth | Melt, acoustic sources, dopants, defects, crystal |
| Transport & Convection (2D) | Convection vs diffusion schematic particle field |

### Panels & Controls

#### Environment
| Control | Effect |
|---------|--------|
| Earth / Space toggle | Switches convection & defect probabilities |
| Temperature (°C) | Affects temperature portion of quality metric |
| Acoustic Intensity | Dampens convection, modulates surface, boosts uniformity |

#### Material
Shows melting point, crystal structure, nominal growth rate, defect probabilities.

#### Growth & Visualization
| Control | Description |
|---------|-------------|
| Growth Progress Slider | Manually sets solid fraction (0–100%) |
| Auto Growth | Time-evolving growth fraction |
| Rotation Speed | Camera spin rate |
| Toggles | Flow lines, dopants, defects, lattice detail |
| Reset | Resets state (growth, defects, dopant field) |

#### Status Panel
Dynamic readouts:
- Crystal Quality (%)
- Defects (count)
- Environment mode
- Dopant Uniformity (%)
- Convection regime
- Acoustic intensity

-

## Dopant Transport Details

1. Grid: Axisymmetric (r,z), size `DOPANT_N × DOPANT_N`.
2. Update per frame:
   - Semi‑Lagrangian advection (Earth only).
   - Explicit diffusion (stability via small dt clamp).
   - Interface band depletion: thickness = `INTERFACE_THICKNESS`.
   - Incorporated dopants sampled angularly to approximate 3D uniformity.

**Uniformity Metric:**
- Angular dispersion (vector resultant length on unit circle).
- Penalized by field standard deviation / mean.
- Acoustic uniformity boost added.

-

## Acoustic Field Heuristics

| Component | Effect |
|-----------|--------|
| Standing-wave superposition | Drives potential signal (visual only) |
| Surface modulation | Low-order harmonic distortion (Y₂₀ + Y₂₂ style) |
| Convection damping | Velocity scaled by `(1 - α·intensity)` |
| Quality boost | Added to uniformity term |

-

## Crystal Quality (Q)

\[
Q = w_u U_{\text{eff}} + w_T T_{\text{score}} + w_A A_{\text{score}} + w_D D_{\text{score}}
\]

Where:
- \(U_{\text{eff}} = \text{uniformity} + \text{acoustic boost}\)
- \(T_{\text{score}}\) penalizes temperature deviation from \(T_m - 50^\circ C\)
- \(A_{\text{score}} = 0.5 + 0.5 I_{acoustic}\)
- \(D_{\text{score}}=1-\text{defect\_probability}\)

All clamped to [0,1].

-

## Defects

| Mode | Distribution Bias |
|------|-------------------|
| Earth | Edge & lower region favored (gravity-related nucleation stress) |
| Space | Uniform random interior |

Regenerated / accumulated depending on version branch. Count displayed in status panel.

-

## Limitations

| Not Modeled | Impact |
|-------------|--------|
| True CFD / Navier–Stokes | Flow patterns idealized |
| Phase-field interface | No morphological instabilities (e.g. dendrites) |
| Thermo-chemical coupling | Temperature not spatially resolved |
| Segregation coefficients (k0) | Incorporation kinetics simplified |
| Stress / Strain | No elastic/plastic defect evolution |
| Real acoustic node solving | Standing wave schematic only |


## Performance Tips

| Action | Benefit |
|--------|---------|
| Reduce `DOPANT_N` (e.g. 48) | Faster PDE updates |
| Disable Dopants/Flow Lines | Fewer GL primitives |
| Lower Acoustic Intensity | Slight reduction in per-frame math |
| Close other GPU apps | Prevent driver contention |


## Troubleshooting

| Symptom | Possible Cause | Fix |
|---------|----------------|-----|
| Crash on tab switch | Material key mismatch | Ensure Unicode saved as UTF‑8 |
| Black 3D viewport | OpenGL context / driver issue | Update drivers / test simple PyOpenGL script |
| No defects appear | Growth still 0% | Increase growth or enable auto growth |
| Uniformity stuck ~100% | Not enough dopant incorporation yet | Allow more growth time |
| Garbled characters | Console encoding | Set `PYTHONUTF8=1` |

-

## Possible Extensions

| Area | Enhancement |
|------|-------------|
| Transport | Full 3D (r,θ,z) or voxel diffusion grid |
| Interface | Phase-field with curvature & anisotropy |
| Thermals | Heat equation + latent heat release |
| Acoustic | Helmholtz or FDTD solver for pressure nodes |
| Defects | Dislocation line objects with motion rules |
| Rendering | Modern OpenGL (VBO/VAO + GLSL shaders) |
| Data Export | CSV / HDF5 snapshots for analysis |

-

## File Layout (Single File Mode)

| Section | Content |
|---------|---------|
| Constants | Physical & visual parameters |
| Materials DB | Quantum material properties |
| Utilities | Clamp, normal calc |
| ConvectionVisualizer | 2D Qt widget (QPainter) |
| CrystalGrowthSimulation | 3D QGLWidget (OpenGL immediate mode) |
| MainWindow | UI & control wiring |
| Entry Point | QApplication bootstrap |

-

## Quick Start Cheat Sheet

| Goal | Action |
|------|--------|
| Run simulation | `python quantumcrystalgrowthearthspaceEnhanced.py` |
| Compare Earth vs Space | Toggle radio buttons & watch dopant flow |
| Improve quality | Set temperature near (Tm − 50°C); raise acoustic intensity; switch to Space |
| See more defects | Use Earth mode + increase growth |
| Reset everything | Click “Reset” |
| Observe uniformity | Watch Dopant Uniformity % while growth proceeds |

-

## License / Attribution

- This simulation’s source code is released under the MIT License.
- Python © Python Software Foundation, under the Python Software Foundation License.
- NumPy © NumPy Developers, under the BSD 3-Clause License.
- Matplotlib © Matplotlib Development Team, under a PSF/BSD-compatible license.
- PyOpenGL © PyOpenGL Developers, under the BSD License.
- pygame © pygame community, under the GNU LGPL 2.1 License.
- OpenGL © Khronos Group, under the OpenGL open standard license.

-

## Support

Open an issue (if in a repo) including:
- OS / Python version
- Console output (traceback)
- Steps to reproduce
- Screenshot (optional)

-

### Disclaimer
This simulation is NOT a substitute for experimental or high-fidelity numerical modeling (CFD, phase-field, FEM). All physical relationships are heuristic approximations for educational visualization.

Enjoy exploring crystal growth physics concepts!
