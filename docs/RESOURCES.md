# Resource log — 4Pi phantom (T1/T2)

**v1.2.0 · 2026-10-30 · Carlos Chacón** (measurements from the v1.0 run)

This log records the memory and time used by the phantom simulation in both views, how much was
left available, and how far the simulation window can be enlarged in the next run, to reduce the
clipping of the field at large angles.

## Machine

| | |
|---|---|
| Computer | Apple Silicon, 8 performance and 4 efficiency cores, **24 GB** RAM |
| Meep | 1.34.0 (conda-forge, mpich); 4 MPI processes per simulation |
| Throughput | 1 process ≈ 43 Mcell·step/s; 4 processes ≈ 98 (2.3×); 6 or 8 processes are slower |
| Measured memory | **139 bytes per Yee cell** (real fields, source included); measured without swapping |
| Free memory at rest | ~10.5 GB (free + inactive) with VS Code and a browser open |

## Run performed (64 illuminations, 32 per view)

| Parameter | Value |
|---|---|
| Object | Registered GT, crop 210 × 266 × 116 voxels (25.5 × 32.3 × 14.1 µm): full support + 6 voxels of medium |
| Resolution | 12 px/µm (step 0.083 µm) |
| Margins | lateral 0.8 µm, walk-off 1.0 µm, z 0.8 µm, PML 1.0 µm |
| Cell | 35–60 Mcells per illumination (mean 49) |
| Memory | max. 8.3 GB estimated; 7.7 GB measured for the 55 M cell (θ = 46°) |
| Time | 0.13 s/step; on average 237 s (bottom) and 273 s (top) per illumination; **4.5 h** in total |

**What did not fit.** At 14 px/µm the 46° cell needed ~14 GB, and with MATLAB open (~7.5 GB) the
machine swapped: the step went from 0.13 to 1.7–2.6 s, 13 to 20 times slower. That is why the
resolution was lowered to 12 px/µm and MATLAB was closed.

**Practical threshold.** Above ~10 GB per simulation the machine swaps with the usual
applications open. The measurements confirm it: 76 M cells took 1.69 s/step and 112 M,
2.65 s/step.

## The observed clipping

The lateral window of the sensor is the object box plus `lateral_margin_um`, extended toward the
illumination side by the walk-off plus `walkoff_margin_um`. Scattered light that leaves through
the sides before reaching the sensor plane enters the PML and is lost. In the 64 simulations, |S|
at the edge of the window reached **3–6 % of its maximum**. That field is lost, and it is the
clipping visible in the reconstructions at large angles.

Angles captured from the base of the object, its deepest part: h = 13.7 µm below the sensor.

- **Side opposite the walk-off:** atan(lateral / h) = **3.3°** with 0.8 µm.
- **Walk-off side:** θ_illumination + atan((walk-off + lateral) / h) = **θ + 7.5°**.

Points closer to the sensor capture larger angles. The loss mainly affects large-angle scattering
from the bottom edge of the object.

## Options for the next run

Computed with the real cell geometry (`python -m fdtd_meep estimate`) on the same 64
illuminations, at 12 px/µm and with `z_margin_um = 0.4`:

| lateral / walk-off (µm) | Cells mean / max (M) | Max memory | Batch of 64 | Opposite side | Walk-off side |
|---|---|---|---|---|---|
| 0.8 / 1.0 (used) | 47.8 / 58.7 | 8.2 GB | 2.1 h | 3.3° | θ + 7.5° |
| 0.8 / 3.0 | 52.8 / 64.1 | 8.9 GB | 2.3 h | 3.3° | θ + 15.5° |
| **1.0 / 4.0** | 56.5 / 68.2 | **9.5 GB** | 2.5 h | 4.2° | **θ + 20.1°** |
| 1.5 / 2.0 | 53.8 / 65.3 | 9.1 GB | 2.3 h | 6.3° | θ + 14.4° |
| 2.0 / 3.0 | 59.1 / 71.1 | 9.9 GB | 2.6 h | 8.3° | θ + 20.1° |
| 1.5 / 6.0 | 64.6 / 77.1 | 10.7 GB | 2.8 h | 6.3° | θ + 28.7° |
| 3.0 / 4.0 | 67.5 / 80.2 | 11.1 GB | 2.9 h | 12.4° | θ + 27.1° |
| 4.0 / 5.0 | 76.4 / 89.8 | 12.5 GB | 3.3 h | 16.3° | θ + 33.3° |

The times assume the machine does not swap. The "batch of 64" does not include the calibration
(2 simulations, ~10 min).

How to read it:

1. **Enlarging the walk-off side is the cheapest option**, because it grows on one side only and
   mostly along one axis. It is also where the clipping was observed, in the large-angle
   illuminations.
2. **Recommended: lateral 1.0 / walk-off 4.0 µm.** It captures almost three times the angle on the
   walk-off side (θ + 20° instead of θ + 7.5°) with 16 % more memory (9.5 GB max.) and about
   25 min more per batch. It fits in the free memory at rest (~10.5 GB) if MATLAB is closed.
3. **With more free memory:** closing VS Code and the browser as well frees another 3–4 GB. Then
   3.0 / 4.0 (11 GB) or 4.0 / 5.0 (12.5 GB) also improve the opposite side.
4. **What does not help:** reducing `z_margin_um` barely changes anything (−2 %), because the
   height of the cell is set by the sensor plane. Going below 12 px/µm leaves fewer than 5 points
   per wavelength and is not validated.
5. **Capturing the whole NA cone (60°) from the base** would require about 23 µm of margin per
   side (~4 times the memory). It is not feasible on this machine without another method, for
   example near-to-far-field surfaces.

Suggested configuration (*Simulation* section of the GUI, or `inputs/config/phantom_4pi.json`):

```json
"simulation": {
  "resolution": 12,
  "lateral_margin_um": 1.0,
  "walkoff_margin_um": 4.0,
  "z_margin_um": 0.4,
  "max_memory_gb": 10.0
}
```

Before launching, check the *Estimate memory* table and the free memory (`memory_pressure` on
macOS). `max_memory_gb` prevents launching an illumination that does not fit.
