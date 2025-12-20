# iDNS — Intelligent Direct Numerical Simulation

Adaptive integration policy for extreme-stiffness dynamical systems via temporal lifting. Verified true DNS (Rε = 1.000) where conventional methods exhibit significant numerical dissipation.

Depending on resolution ratio R = N/(1.6√Re), the method operates as:
- **Spectral DNS** (R ≥ 1): All dynamically active scales resolved; dissipation purely viscous
- **DNS-Coarse** (R < 1): Exact Galerkin truncation; large and medium-scale structures captured with precision while sub-Kolmogorov scales are truncated by spectral projection rather than artificial dissipation

Unlike LES or SGS methods, DNS-Coarse introduces no modeled dissipation — retained modes (k ≪ kmax) evolve under exact Navier-Stokes dynamics. Energy exits at the truncation boundary, not through numerical damping. Adaptive integration framework for extreme-stiffness dynamical systems via temporal lifting. Verified true DNS (Rε = 1.000) where conventional methods exhibit significant numerical dissipation. At under-resolved grids (R < 1), operates as DNS-Coarse: exact Galerkin truncation preserving large and medium-scale structures while sub-Kolmogorov modes exit via spectral projection — no artificial dissipation, no subgrid models. Validated against NASA Glenn WRLES (DeBonis, 2013) and spectral DNS benchmarks (Chandler & Kerswell, 2013).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17730872.svg)](https://doi.org/10.5281/zenodo.17730872)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
Mathematics proof and LEAN verification uploaded for temporal lifting, see zip file "NavierStokes-GlobalRegularity-Lean4.zip" and includes readme how to verify. 

**Version:** 2.0.0  
**Author:** Jeffrey Camlin  
**ORCID:** [0000-0002-5740-4204](https://orcid.org/0000-0002-5740-4204)  
**Date:** 14 December 2025  
**Organization:** Red Dawn Academic Press & AI Lab, Milwaukee WI

<p align="center">
  <img src="https://i0.wp.com/reddawnacademicpress.org/wp-content/uploads/2025/11/teaser-63.png?w=913&ssl=1" alt="iDNS Taylor-Green Vortex" width="700">
</p>

## Repository Structure
```
iDNS/
├── LICENSE
├── README.md
├── 1_solvers/
│   ├── idns_v10_taylor_green_public.py    # 3D Taylor-Green vortex
│   └── idns_v10_kolmogorov_public.py      # 2D Kolmogorov flow
├── 2_utilities/
│   ├── idns_visualizer_v2.py              # Snapshot visualization
│   ├── extract_spectrum.py                 # Energy spectrum extraction
│   └── npz_inspect.py                      # Snapshot inspection
├── 3_data/
│   ├── taylor-green-10-5-data.zip
│   ├── taylor-green-10-8-data.zip
│   ├── idns-nasa-debonis_validation.zip
│   ├── kolmogorov-validation-chandler-kerswell.zip
│   ├── kolmogorov-invariance.zip
│   └── idns-Kolmogorov Flow-xtreme-N512_Re100000000.zip
└── 4_charts-and-visuals/
    ├── TGV_Re10-5_diagnostics.pdf
    ├── TGV_Re10-8_4panel_XZ_t5-N512.png
    ├── TGV_Re10-8_N512_3planes_t5.png
    ├── bkm-taylor-green-10-5-and10-8.pdf
    ├── energy-spectra-taylor-green-10-5-and-10-8.pdf
    └── ...
```

## idns_v10_taylor_green.py

3D Taylor-Green vortex solver with temporal lifting. Validated against DeBonis (2013) NASA Glenn WRLES benchmark. Achieves Rε = 1.000 (true DNS) where benchmark reports Rε = 1.46.

**Hardware Requirements:**
- N = 128³: Runs on laptop/desktop with 8 GB RAM
- N = 256³ or higher: Requires additional RAM (16+ GB recommended)

**DNS vs DNS-Coarse:**
The resolution ratio R = N/(1.6√Re) determines the operating regime:
- **R ≥ 1 (Spectral DNS):** All scales resolved; dissipation purely viscous
- **R < 1 (DNS-Coarse):** Exact Galerkin truncation; large/medium structures preserved, sub-Kolmogorov modes truncated — no artificial dissipation

For Taylor-Green, true DNS threshold is approximately N = 512³ at Re = 10⁵, (see idns-resolution-ratio-chart.pdf in charts-and-visuals)

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--N` | `128` | Grid resolution (N³) |
| `--Re` | `100000` | Reynolds number |
| `--t_end` | `5.0` | Final simulation time |
| `--dt_base` | `1e-4` | Base computational timestep |
| `--snapshots` | off | Save velocity field snapshots (.npz, large files) |
| `--outdir` | `results` | Output directory |
| `--print_every` | `100` | Console output interval |

### DeBonis Benchmark Validation (Re = 1600, t = 0 → 20)
```bash
python idns_v10_taylor_green.py --N 64 --Re 1600 --t_end 20.0
python idns_v10_taylor_green.py --N 128 --Re 1600 --t_end 20.0
python idns_v10_taylor_green.py --N 256 --Re 1600 --t_end 20.0
python idns_v10_taylor_green.py --N 512 --Re 1600 --t_end 20.0
```

All resolutions achieve R_ε = 1.000 ± 0.002, confirming zero numerical dissipation. The 256³ run completes in ~7 days on a laptop; 512³ requires cloud compute (~2 hr on 64-core Azure).

### High-Reynolds Cascade (Re = 10⁵)
```bash
python idns_v10_taylor_green.py --N 128 --Re 1e5 --t_end 5.0
```

Peak vorticity amplifies 14× through the cascade. BKM integral remains bounded (37.1). Runtime: ~80 hours on dual-core laptop.

### Extreme Stiffness (Re = 10⁸)
```bash
python idns_v10_taylor_green.py --N 256 --Re 1e8 --t_end 5.0
```

DNS-Coarse regime (R = 0.016). Captures exact Galerkin dynamics with large-scale fidelity. Funnel-shaped vortex cores visible during sheet-to-tube rollup.

### With Velocity Snapshots
```bash
python idns_v10_taylor_green.py --N 128 --Re 1e5 --t_end 5.0 --snapshots
```

Saves .npz files at t = 3.0, 4.0, 4.5, 5.0 containing full velocity fields (ux, uy, uz). Warning: large files (~500 MB each at 128³).

---

## idns_v10_kolmogorov.py

2D Kolmogorov flow solver with temporal lifting and optional Cascade-Completion Clamping (CCC). Validated against Chandler & Kerswell (2013) spectral DNS benchmark.

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--mode` | `idns` | Integration mode: `conservative` (fixed dt) or `idns` (temporal lifting) |
| `--Re` | `20000` | Reynolds number |
| `--N` | `256` | Grid resolution (N×N) |
| `--t_end` | `100.0` | Final simulation time |
| `--seed` | `48` | Random seed for initial perturbation |
| `--ccc` | `off` | Cascade-Completion Clamping: `off`, `periodic`, or `adaptive` |
| `--ccc_interval` | `200` | Steps between CCC (periodic mode) |
| `--ccc_threshold` | `2.0` | Pile-up factor threshold (adaptive mode) |
| `--outdir` | `results` | Output directory |

### Example: Re = 10⁸
```bash
python idns_v10_kolmogorov.py --mode idns --Re 1e8 --N 512 --t_end 50.0
```

Temporal lifting alone stabilizes integration at extreme Reynolds numbers without artificial dissipation.

### Reynolds Invariance (Re = 2×10⁴ vs Re = 10⁶)
```bash
python idns_v10_kolmogorov.py --mode idns --Re 2e4 --N 256 --t_end 100.0
python idns_v10_kolmogorov.py --mode idns --Re 1e6 --N 256 --t_end 100.0
```

Despite 50× viscosity separation, energy and enstrophy statistics match within 3%, confirming Reynolds-invariant controller response. Dissipation scales exactly with viscosity: ε₂₀ₖ/ε₁ₘ = 48.7 vs theoretical 50.0.

### Chandler & Kerswell Validation (F=1, kf=4, N=256, t=1000)
```bash
python idns_v10_kolmogorov.py --mode idns --Re 40 --N 256 --t_end 1000.0
python idns_v10_kolmogorov.py --mode idns --Re 60 --N 256 --t_end 1000.0
python idns_v10_kolmogorov.py --mode idns --Re 80 --N 256 --t_end 1000.0
python idns_v10_kolmogorov.py --mode idns --Re 100 --N 256 --t_end 1000.0
```

All cases achieve R_ε = 1.0 ± 0.03, confirming exact energy budget closure with zero numerical dissipation. Resolution ratio R ≫ 1 places these firmly in the DNS regime.

---

## idns_visualizer.py

Visualization utilities for Taylor-Green vortex snapshots. Generates 2D slices and optional 3D isosurface renders.

### Snapshot Files

The visualizer requires `.npz` snapshot files containing velocity field data (~48 MB each at 128³). Due to size constraints, snapshots are not included in the repository.

To generate snapshots, run the Taylor-Green solver with the `--snapshots` flag:
```bash
python 1_solvers/idns_v10_taylor_green_public.py --N 128 --Re 1e5 --t_end 5.0 --snapshots
```

This outputs `field_t3.0.npz`, `field_t4.0.npz`, `field_t4.5.npz`, `field_t5.0.npz` in the `results/` directory, ready for visualization.

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `npz_path` | (required) | Path to .npz snapshot file |
| `--Re` | `10^5` | Reynolds number label for plot titles |
| `--iso` | off | Generate 3D isosurface (slower) |
| `--iso_level` | 25% of max | Isosurface vorticity threshold |

### Basic Usage
```bash
python idns_visualizer.py field_t4.5.npz
```

### With Reynolds Label
```bash
python idns_visualizer.py field_t4.5.npz --Re 1e5
python idns_visualizer.py field_t4.5.npz --Re 1e8
```

### With 3D Isosurface
```bash
python idns_visualizer.py field_t4.5.npz --Re 1e8 --iso
```

### Output Files
```
1-TG_composite_4panel.png    — XZ slice: velocity, vorticity, composites
2-TG_funnel_detail.png       — XZ vorticity with/without contours
3-TG_composite.png           — 3-panel velocity + vorticity contours
4-TG_velocity_slices.png     — XY, XZ, YZ velocity magnitude
5-TG_vorticity_slices.png    — 6 z-slices of vorticity
6-TG_vorticity_3panel.png    — XY, XZ, YZ vorticity magnitude
7-TG_isosurface_colored.png  — 3D vortex tubes (with --iso)
```

---

## npz_inspect.py

Quick inspection utility for iDNS snapshot files. Displays keys, shapes, dtypes, and value ranges.
```bash
python npz_inspect.py field_t4.5.npz
```

Output:
```
File: field_t4.5.npz
Keys: ['ux', 'uy', 'uz', 't']
  ux: shape=(128, 128, 128), dtype=float64, min=-1.234, max=1.456
  uy: shape=(128, 128, 128), dtype=float64, min=-1.198, max=1.203
  uz: shape=(128, 128, 128), dtype=float64, min=-1.301, max=1.287
  t: 4.5001234

npz_inspect v10.0.0 complete.
```

---

## Repository Contents
```
iDNS/
├── README.md
├── idns_v10_taylor_green.py    — 3D Taylor-Green solver
├── idns_v10_kolmogorov.py      — 2D Kolmogorov solver
├── idns_visualizer.py          — Snapshot visualization
├── npz_inspect.py              — Snapshot inspection utility
├── data/
│   └── field_t*.npz            — Pre-generated snapshots
└── results/                    — Output directory (created on run)
```

---

## Citation
```bibtex
@software{camlin2025idns,
  author = {Camlin, Jeffrey},
  title = {iDNS: Intelligent Direct Numerical Simulation},
  year = {2025},
  doi = {10.5281/zenodo.17730873},
  url = {https://github.com/RED-DAWN-AI/iDNS}
}
```

---

## References

### Primary Sources

Camlin, J. (2025). iDNS: A state-dependent temporal policy for extreme-stiffness dynamical systems, with application to the 3-D Taylor–Green vortex at Re = 10⁵. *arXiv preprint* (forthcoming). Version 2 archived on PhilArchive: [CAMIIA-3.pdf](https://philpapers.org/archive/CAMIIA-3.pdf)

Camlin, J. (2025). Global regularity for Navier–Stokes on T³ via bounded vorticity–response functionals. *The Scholarly Journal of Post-Biological Epistemics*, 1(2), 1–14. [doi:10.63968/post-bio-ai-epistemics.v1n2.012](https://doi.org/10.63968/post-bio-ai-epistemics.v1n2.012)

Camlin, J. (2025). Temporal lifting as latent-space regularization for continuous-time flow models in AI systems. *arXiv:2510.09805* [cs.LG]. [doi:10.48550/arXiv.2510.09805](https://doi.org/10.48550/arXiv.2510.09805)

### Benchmarks & Background

DeBonis, J. R. (2013). Solutions of the Taylor-Green vortex problem using high-resolution explicit finite difference methods. NASA/TM-2013-217850.

Chandler, G. J. & Kerswell, R. R. (2013). Invariant recurrent solutions embedded in a turbulent two-dimensional Kolmogorov flow. *J. Fluid Mech.* 722:554–595.

van Rees, W. M., Leonard, A., & Pullin, D. I. (2011). A comparison of vortex and pseudo-spectral methods. *J. Comput. Phys.* 230(8):2794–2805.

Beale, J. T., Kato, T., & Majda, A. (1984). Remarks on the breakdown of smooth solutions for the 3-D Euler equations. *Commun. Math. Phys.* 94:61–66.
