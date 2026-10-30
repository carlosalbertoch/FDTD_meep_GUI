# FDTD_meep — Technical documentation

**v1.2.0 · 2026-10-30 · Carlos Chacón**

> v1.2: graphical interface ([GUIDE.md](GUIDE.md)), free mode (`free`), illumination lists,
> sphere objects, medium index taken from the object (no index correction), and automatic crop.
> The solver and the v1.1 steps are unchanged.

1. [What it does](#what-it-does)
2. [How it works](#how-it-works)
3. [Conventions](#conventions)
4. [Input parameters](#input-parameters)
5. [Commands and outputs](#commands-and-outputs)
6. [Your own objects from MATLAB](#your-own-objects-from-matlab)
7. [Validation](#validation)
8. [Known limitations](#known-limitations)
9. [Modules](#modules)

---

## What it does

It reproduces a transmission TDM acquisition with FDTD. For each selected illumination of a real
acquisition, it computes the field scattered by the object, carries it to the camera through the
microscope optics, and forms the off-axis hologram. The result can go into the reconstruction
pipeline exactly like the experimental holograms. Several views are supported: in a 4Pi system,
the second view is equivalent to the sample rotated by 180°. The backscattered light recorded by
the opposite camera (4Pi reflection) is also simulated.

Without real data, the free mode simulates any object for a list of illumination directions and
saves the field at the sensor.

## How it works

```text
GT / object n(x,y,z) ──prepare──► object.h5
real holograms ───────select───► specular peaks (k_RAW) + FPS or list selection + real_subset/
                     calibrate ─► object → camera transform (per view)
                     run ───────► Meep (scattered field) → pupil → camera → hologram → sim_subset/
                     reflect ───► backscatter → opposite camera → reflection holograms
                     compare ───► demodulated real vs simulated fields
object + list ────────free ────► Meep (scattered field) → fields at the sensor (no camera)
```

### 1. Object (`prepare`)

- `source = "registered_gt"`: the registered ground truth (MATLAB v7.3 file with the variable
  `phantom_gt_dn` in V(y,x,z) convention, and `meta.rec_voxel_um_xyz`) is cropped to `crop_xyz`.
  The index is n = n_medium + Δn, **with no correction**. n_medium is the one stored in the file
  (`meta.n_medium`, or `meta.phantom_n_immersion` if that is missing) unless `object.n_medium` is
  set. Without `crop_xyz`, the crop is automatic: the bounding box of the support
  (|Δn| > `meta.mask_threshold`, 1e-4 by default) plus `crop_pad_um` on every side. The cropped
  volume is centred at (0, 0, 0) in the simulation; `calibrate` gives the alignment with each
  camera. (v1.1 added (n_immersion_gt − n_medium) with a soft mask; this is applied only if
  `n_immersion_gt` is set explicitly.)
- `source = "voxel_h5"`: your own volume in the exchange format (see
  [Your own objects](#your-own-objects-from-matlab)).
- `source = "sphere"`: a sphere with a smooth tanh edge (`radius_um`, `dn`, `voxel_um`,
  `edge_um`).

`prepare` prints the n_medium in use and where it comes from, warns if it differs from N0 of
`config_manip.txt`, and warns if the object touches the border of the crop.

### 2. Specular peaks and selection (`select`)

Each iNNN.pgm hologram (`Intensite_ref.pgm` is excluded) is processed as follows: Tukey window →
FFT → patch of the +1 order, of size dimHolo = 2·ceil(R_Ewald·NA/N0), centred according to
CIRCLE_CX/CY → maximum → subpixel centroid. The position of the peak relative to the patch centre
gives the illumination vector in RAW convention.

By default, the illuminations are ordered by **farthest-point sampling** in (kx, ky), starting
with the most central one, and the first `n` are taken. These `n` are always a uniform sample of
the scan pattern. With `selection.lists`, the order and the holograms come from a list file
instead (see [Illumination lists](#illumination-lists)). The selected holograms are copied,
renumbered (i000…), to `real_subset/`.

### 3. Scattered-field formulation (core, `scattered.py`)

Subtracting the Maxwell equations of the incident field from those of the total field gives

```
∂(ε E_s)/∂t = ∇×H_s − J,     J = ∂P/∂t,     P = (ε − ε_medium)·E_inc
```

so the scattered field E_s is obtained without simulating the incident wave:

- **Discrete incident wave.** k_t is the continuous value measured by the camera; κ_z comes from
  the Yee dispersion relation, and the polarization is perpendicular to the discrete k̃. The
  source is therefore exact on the grid.
- **Source.** For each component (Ex, Ey, Ez), Δε is evaluated at its Yee position with the same
  linear interpolation in ε that `MaterialGrid` uses, and multiplied by −iω·e_c·E_inc.
- **Gaussian pulse and DFT** at the centre frequency, normalized by the actual spectrum of the
  applied current (Σ J(t)e^{iωt}Δt/√2π). Real fields are used, with half the memory; the result
  is identical to complex fields.
- **Cell.** It contains the object box plus `lateral_margin_um`, plus `z_margin_um`, plus the
  sensor plane and the PML. With oblique illumination, the transmitted field leaves the object
  shifted by (z_sensor − z_min)·tan θ, so the cell is extended **only toward that side** by that
  shift plus `walkoff_margin_um`. The number of cells per axis is even, to keep the Yee grid
  aligned with the origin.
- **Stop.** The run ends when |E|² at the sensor has dropped by `field_decay_by` from its maximum,
  once the source is off, or when `max_runtime_periods` is reached.
- **Calibration.** An empirical modulus factor, measured against Mie and dependent on the
  resolution, is applied (see [Validation](#validation)).

Each illumination runs in a new process (`mpirun -np mpi_processes`), so memory is always
released.

### 4. Optics and camera

S = E_s/E_inc (relative scattered field) is multiplied by the illumination ramp and passed
through the NA pupil. It is then propagated with the angular spectrum to the focused plane, and
the incident wave is propagated too. Finally it is interpolated onto the camera grid
(dx = TPCAM/Gtot) with the calibrated transform r_cam = M·r_obj + t.

The hologram is I = |U_cam + U_ref|². The carrier is chosen so that the U·U_ref* term (or
U*·U_ref if the view uses the −1 order) falls at the centre of the pipeline patch. The mean level
is matched to that of the real hologram and quantized to 8 bits.

### 5. Calibration with real data (`calibrate`)

1. Two normal incidences are simulated: light entering through the **zmin** face (object as is)
   and through the **zmax** face (object rotated 180° about y).
2. For each view, the most central real hologram is demodulated: patch, recentred peak and
   pupil. The phase and log-amplitude background is removed with a 2D polynomial fitted on an
   outer annulus, because it contains the beam curvature and aberrations.
3. The best match with the simulated field is searched for. The 8 symmetries of the square (M),
   U or U* (order +1 or −1), the shift by cross-correlation (t) and the focal plane, by a z
   scan, are all tried.
4. With two views, they must enter through opposite faces. The pairing with the highest sum of
   complex correlation plus amplitude correlation is chosen.

You never flip the object by hand for the second view: this step decides which face each view
uses.

### 6. Reflection (`reflect`)

In the 4Pi each view has its own reference arm, so backscattered light also forms holograms. The
reflection recorded by the top camera uses the illuminations of the *bottom* transmission (the
source is on top), and the one recorded by the bottom camera uses those of the *top*
transmission:

| Reflection view | Illuminations (selection and calibration) | Camera (calibration) |
|---|---|---|
| `refl_top` | `bottom` | `top` |
| `refl_bottom` | `top` | `bottom` |

1. **Simulation.** As in transmission, with the object in the orientation of the illumination
   view and the same k_obj, but the sensor plane is `sensor_gap_um` **below** the object. With
   the scattered-field formulation there is no incident wave at the sensor: only the backscatter
   is recorded. The walk-off has the same sign as in transmission, since reflection preserves
   the transverse k. The PML on the +z side, where the much stronger forward field leaves, can be
   made thicker with `pml_far_um`. 1 µm is enough: 2 µm does not change the result.
2. **Transfer to the opposite camera.** The camera was calibrated with the object rotated 180°
   about y relative to the illumination view. In that frame, x → −x and z → −z, the backscatter
   travels toward +z, and the pupil, focus, M, t and ±1 order of that camera are applied
   unchanged. The transverse k at the camera is k_c = M·(−kx, ky).
3. **Holograms.** U = S·exp(i k_c·r), **without the 1 of the transmitted wave**, and
   I = |U + r_p·exp(i k_c·r) + R|². r_p is an optional parasitic specular, with amplitude relative
   to the incident one; the true specular, if present, falls at k_RAW = s·k_c/k. |R| is the same
   for the whole view: `reference_ratio` × the 99.9th percentile of |U|, median over the view.
   The gray scale is also one per view, so that the mean level matches the real holograms of
   that camera.

Outputs in `outputs/<run>/<reflection_view>/`:

- `fields/iNNN.npz`: fields at the sensor, in the same format as in transmission, plus
  `reflection`.
- `camera_fields/iNNN.npz`: `spectrum_patch`, the complex +1-order patch of the cross term U·R*,
  normalized by |R|, with no noise or quantization. It includes the expected specular position
  (`drow_specular`, `dcol_specular`, in bins from the patch centre), `fmaxHolo` and `REwald`. It
  is the reference for testing Rytov normalizations without a specular.
- `sim_subset/`: holograms without a parasitic specular; `sim_subset_spec<r>/`: with r_p = r.
  Each folder contains the camera's `config_manip.txt` and `recon.txt`, `Intensite_ref.pgm`,
  `seleccion.txt`, `sim_report.csv` and `sim_parameters.json`. For each illumination,
  `sim_report.csv` gives the expected specular position, the patch maximum and its distance to
  the specular (`peak_to_specular_bins`), the RMS and maximum of |U|, and saturation.
- `<view>_resumen.png`: hologram, +1-order spectrum (with the expected specular marked), and
  amplitude and phase of U, for θ ≈ 0°, 12°, 25° and 40°.

### 7. Lateral reference figures (`reference`)

An optional feature, separate from the campaign: nothing changes if it is not run. It repeats the
**central illumination** (kx = ky = 0) with DFT monitors on XZ and/or YZ slices of the whole
domain without PML. It shows the field **inside** the FDTD domain, not derived from the sensor or
from `camera_fields`.

- A transmission view and the reflection view with the same illumination (`bottom` and
  `refl_top`; `top` and `refl_bottom`) share the same simulation; only the camera changes. Two
  simulations (~6 GB, ~7 min each) cover the four cases.
- `z_pad_um` adds margin in z on both sides, to see the wave before and after the object.
- The figures are drawn in the laboratory frame: +z is the top of the microscope in every view.
  The arrow marks k_inc and the label, the camera that records.
- Total field \|E_s + E_inc\|² and scattered field \|E_s\|², summing the three components and
  normalizing by \|E_inc\|². With `compressed_arcsinh`, the image is arcsinh(αI)/arcsinh(αI_max).
  With automatic α, α = 1/(typical I on the camera side, outside the object): linear up to that
  level and logarithmic above it. In reflection, α ≈ 10⁴ shows the backscatter next to the much
  stronger forward field.
- The contour is the half-height isoline of Δn of the simulated object.

```bash
./.env/bin/python -m fdtd_meep reference --config $CFG                  # all 4 views
./.env/bin/python -m fdtd_meep reference --config $CFG --view refl_top  # one view
```

Outputs: `outputs/<run>/reference_views/fdtd_<orientation>_incidencia_normal.npz`, with the raw
slices (`slice_<plane>_{x,y,z}`, `slice_<plane>_Es_*`, `slice_<plane>_Einc_*`, simulation frame);
`<view>/reference_views/lateral_field_{total,scattered}_<plane>_<mode>_<cmap>.{png,svg}`; and
`reference_views/lateral_reference_4casos_*.{png,svg}` (2 × 2: transmission and reflection ×
top and bottom). In the SVG the map is rasterized; axes, text and colorbar remain editable.

### 8. Free mode (`free`)

Without real holograms or calibration. The object (sphere, `voxel_h5` or `registered_gt`) is
simulated as is (orientation zmin, light entering through the z-min face) for every direction of
the list in `illuminations`, in transmission and/or reflection (`modes`). (kx, ky) are in the
object frame. The output is the field at the sensor, in the same `.npz` format as `run`/`reflect`,
with a report and a summary figure (|1 + S| or |S| and phase at the sensor).

## Conventions

| Aspect | Convention |
|---|---|
| Units | µm; λ in vacuum; k in rad/µm |
| Time | exp(−iωt) (Meep): a forward-travelling wave is exp(+ik·r) |
| Volumes | Python `n[z, y, x]`; MATLAB `V(y, x, z)`, as the registered GT |
| 2D arrays | `[row = y, column = x]` |
| Experimental k | RAW: +dcol → +kx, +drow → +ky (camera of each view); kx_norm = kx/(2π·n_medium/λ) |
| k in Meep | k_obj = Mᵀ·(s·k_RAW), with s = −1 if the view uses the −1 order |
| Orientations | `zmin`: light enters through the minimum-z face of the object; `zmax`: object rotated 180° about y |
| Fourier patch | dimHolo = 2·ceil(R_Ewald·NA/N0); its centre is 1 bin after CIRCLE_CX/CY |

## Input parameters

A single JSON file (`inputs/config/*.json`). Relative paths are resolved from the repository
root. An unknown parameter raises an error. The GUI writes this file.

### General

| Key | Description |
|---|---|
| `name` | Name of the run |
| `output` | Results folder (default `outputs/<name>`) |
| `views` | `{"view": "folder"}`. Each folder contains `config_manip.txt`, `iNNN.pgm` and `Intensite_ref.pgm`. NA, N0, LAMBDA, TPCAM, F_TUBE, F_OBJ, RF, DIM_ROI and CIRCLE_CX/CY are read from `config_manip.txt`. Without `views`, the configuration is a free-mode one |
| `wavelength_um` | Free mode only: λ in vacuum (µm) |
| `illuminations` | Free mode only: illumination list (CSV, see [Illumination lists](#illumination-lists)) |
| `modes` | Free mode only: `["transmission"]`, `["reflection"]` or both |

### `object`

| Key | Default | Description |
|---|---|---|
| `source` | `registered_gt` | `registered_gt`, `voxel_h5` or `sphere` |
| `file` | `inputs/phantom/phantom_ground_truth_registered.mat` | Registered GT or `.h5` volume |
| `crop_xyz` | automatic | `[[x0,x1],[y0,y1],[z0,z1]]`: the faces of the box in GT voxels, 0-based with the end excluded (MATLAB: `V(y0+1:y1, x0+1:x1, z0+1:z1)`). Must contain the whole support plus some medium |
| `crop_pad_um` | 0.8 | Medium added around the support in the automatic crop (µm per face) |
| `n_medium` | the object's | Medium index. By default, the one stored in the file (`.mat` `meta.n_medium` / `meta.phantom_n_immersion`; `.h5` attribute `n_medium`); if the object does not store it, N0 of `config_manip.txt`. `prepare` warns if it differs from N0 |
| `n_immersion_gt` | null (no correction) | v1.1 compatibility only: Δn referred to another immersion |
| `soft_mask_dn` | 0.01 | Scale of the soft mask of that correction |
| `center_um` | `[0,0,0]` | Centre of the object in the simulation |
| `radius_um`, `dn`, `voxel_um`, `edge_um` | 1.5, 0.03, 0.1, 0.05 | `sphere` only: radius, Δn, voxel and width of the smooth tanh edge |

### `simulation`

| Key | Default | Description |
|---|---|---|
| `resolution` | 12 | px/µm. Calibrated at 12, 14 and 20. Below 12 there are fewer than 5 points per wavelength in n = 1.55 |
| `pml_um` | 1.0 | PML thickness |
| `pml_far_um` | null | Reflection only: PML thickness on +z, the side of the forward field (null = `pml_um`) |
| `lateral_margin_um` | 0.8 | Lateral margin around the object (x, y) |
| `walkoff_margin_um` | 1.0 | Additional margin toward the walk-off side |
| `z_margin_um` | 0.4 | Margin in z below and above the object box |
| `sensor_gap_um` | 0.4 | Distance from the sensor plane to the object box |
| `pulse_fwidth` | 0.4 | Relative spectral width of the pulse (does not change the result; smaller = slower) |
| `field_decay_by` | 1e-8 | Decay of \|E\|² at the sensor that stops the simulation |
| `max_runtime_periods` | 500 | Safety limit, in periods |
| `polarization` | `Ex` | `Ex` or `Ey` in the object frame (projected ⟂ k) |
| `mpi_processes` | 4 | Processes per simulation (1 = no mpirun). On an M-series Mac with 8 performance cores, 4 is optimal |
| `max_memory_gb` | null | If the estimated memory exceeds it, the illumination is not launched |
| `bytes_per_cell` | 139 | Memory per Yee cell used in the estimates (measured; see [RESOURCES.md](RESOURCES.md)) |

### `selection`

| Key | Default | Description |
|---|---|---|
| `n` | 32 | Illuminations per view |
| `tukey` | 0.10 | Tukey window of the extractor (the same as the pipeline) |
| `lists` | `{}` | `{"view": "list.csv"}`: the selection follows the list instead of FPS (`n` is ignored) |

### Illumination lists

A CSV file with a header row; `#` starts a comment. One row per illumination, in simulation order,
with **one** of these column sets: `kx,ky` (direction cosines kx/k, ky/k, with
k = 2π·n_medium/λ and kx² + ky² < 1), `theta_deg,phi_deg` (kx = sin θ cos φ, ky = sin θ sin φ) or
`filename` (a real hologram, 4Pi only). `name` is optional (default i000, i001…). The
`illumination_positions.csv` written by `select` is also accepted.

- **Free mode**: (kx, ky) in the object frame; light enters through the zmin face.
- **4Pi** (`selection.lists`): (kx, ky) in the RAW convention of the view's camera. Each
  direction is assigned to the real hologram with the closest specular peak, with no repetitions;
  with `filename`, that hologram is used. The selection is written exactly as with FPS
  (`seleccion.txt`, `real_subset/`…), so `calibrate`, `run` and `reflect` are unchanged.

### `calibration`

| Key | Default | Description |
|---|---|---|
| `z_range_um` | `[-14, 14]` | Search range of the focal plane (simulation frame) |
| `z_step_um`, `z_refine_um` | 1.0, 0.25 | Coarse and fine steps of that search |
| `background_annulus_um` | `[17, 33]` | Annulus (radius in µm from the FOV centre) where the background is fitted; must exclude the object |
| `background_poly_order` | 4 | Order of the background polynomial |
| `match_radius_um` | 18 | Radius of the amplitude correlation |

### `reflection`

| Key | Default | Description |
|---|---|---|
| `views` | `refl_top` ← illumination `bottom`, camera `top`; `refl_bottom` ← illumination `top`, camera `bottom` | Reflection views. Illumination and camera must enter through opposite faces |
| `reference_ratio` | 2.0 | \|R\| / level of \|U\| (99.9th percentile, median over the view) |
| `parasitic_specular` | `[]` | Amplitudes r_p of the parasitic specular; each one produces `sim_subset_spec<r>/` |

### `reference_views`

| Key | Default | Description |
|---|---|---|
| `enabled` | false | If true, `run` and `reflect` also produce the figures of their view. The `reference` command always works |
| `planes` | `["YZ"]` | `XZ`, `YZ` or both |
| `position_um` | null | Slice position (x for YZ, y for XZ) in the laboratory frame; null = object centre |
| `field_modes` | `["total","scattered"]` | Fields to draw |
| `intensity_mode` | `compressed_arcsinh` | or `raw` |
| `arcsinh_alpha` | null | Fixed α; null = automatic |
| `colormap` | `turbo` | or `viridis` |
| `overlay_contour` | true | Phantom contour |
| `save_png`, `save_svg` | true, true | Formats |
| `z_pad_um` | 3.0 | Extra margin in z on each side (this simulation only) |

## Commands and outputs

`python -m fdtd_meep <command> --config <json>`. Outputs in `outputs/<run>/`:

| Command | Writes |
|---|---|
| `prepare` | `object/object.h5` |
| `select` | `<view>/selection/` (`illumination_positions_all.csv`, `illumination_positions.{csv,mat}`, `seleccion.txt`, `seleccion_k.png`) and `<view>/real_subset/` |
| `estimate` | A table on screen: cells and memory per illumination (uses the calibration if it exists) |
| `calibrate` | `calibration/center_{zmin,zmax}.npz`, `calibration.json`, `calibration.png` |
| `run --view V [--limit N]` | `<view>/fields/iNNN.{npz,log}` (fields at the sensor) and `<view>/sim_subset/` (iNNN.pgm, `Intensite_ref.pgm`, `config_manip.txt`, `recon.txt`, `seleccion.txt`, `sim_report.csv`, `sim_parameters.json`). Resumable: existing `.npz` files are skipped |
| `compare` | `comparison/comparison.{png,json}` |
| `reflect --view R [--limit N]` | `R/fields/`, `R/camera_fields/`, `R/sim_subset*/`, `R/R_resumen.png` (see [Reflection](#6-reflection-reflect)). Resumable |
| `free [--mode M] [--limit N]` | Free mode, no real holograms: `M/fields/<name>.npz` (same format), `M/free_report.csv`, `M/free_parameters.json`, `M/illuminations_used.csv`, `M/M_resumen.png` (\|1 + S\| or \|S\| and phase at the sensor). M = `transmission` or `reflection`. Resumable |
| `reference [--view V]` | Lateral reference figures (see §7) |
| `example-sphere` | `outputs/example_sphere/` |

Contents of each `fields/iNNN.npz`:

- `Es_{Ex,Ey,Ez}` and `Einc_{Ex,Ey,Ez}`: fields at the sensor, `[y, x]`;
- `x`, `y`, `z`: coordinates;
- `kx`, `ky`, `kappa`: wave vector;
- `pol`: polarization;
- `cell`, `cell_center`, `walkoff`, `n_cells`: cell geometry;
- `meep_time`, `resolution`, `orientation`.

`sim_report.csv` gives Δp for each illumination: the distance in bins between the simulated and
the real specular peak.

Every output file is described in [OUTPUTS.md](OUTPUTS.md).

## Your own objects from MATLAB

`matlab/exportar_objeto_meep.m` writes the exchange format from a volume in the GT convention,
V(y, x, z):

```matlab
VOX = 0.1;                                   % µm
g = ((0:59) - 29.5) * VOX;
[X, Y, Z] = meshgrid(g, g, g);               % arrays (y, x, z)
dn = 0.02 * (X.^2/2^2 + Y.^2/1.5^2 + Z.^2/1^2 < 1);
exportar_objeto_meep(dn, VOX, 1.500, 'inputs/objects/ellipsoid.h5');   % Δn → n = 1.5 + Δn
```

Format: `/n` is the absolute index, which h5py reads as `n[z, y, x]` (MATLAB writes it as
n(x, y, z)); `/coords/{x,y,z}_um` contain the voxel centres, centred at 0; the `axis_order`
attribute is `'zyx'`, and the `n_medium` attribute is the medium index that the simulation uses.
To use it, in the configuration:
`"object": {"source": "voxel_h5", "file": "inputs/objects/ellipsoid.h5"}`.

To set the medium index of the registered GT (`.mat`), add `meta.n_medium`; `phantom_gt_dn`
stays Δn relative to that medium:

```matlab
load(f, 'meta');  meta.n_medium = 1.5;  save(f, 'meta', '-append');
```

Recommendations:

- smooth edges, for example 0.5·(1 − tanh(d/w)) with w ≈ 0.05 µm;
- the border of the volume should equal n_medium (the function warns otherwise);
- voxel ≤ 0.12 µm.

## Validation

- **Scattered field vs exact Mie.** The 2D complex angular spectrum of the field scattered by a
  sphere (Δn = 0.03, R = 1.5 µm) is compared at three angles up to 35°. After the modulus
  calibration, the complex residual is 0.9–1.6 % at 20 px/µm, 1.6–2.6 % at 14 px/µm and
  2.2–3.7 % at 12 px/µm, with a phase ≤ 0.06 rad. Modulus calibration (Mie/FDTD): 1.070 at
  20 px/µm, 1.056 at 14 px/µm and 1.19 ± 3 % at 12 px/µm. Its origin has not been identified; it
  is applied as an empirical correction.
- **Real vs simulated fields (4Pi phantom, v1.1 run).** With 32 illuminations per view, the
  simulated specular peak falls on average 0.12 bins from the real one (maximum 0.25). Complex
  correlation of the demodulated fields:

  | View | 0.7° | 12° | 41° |
  |---|---|---|---|
  | bottom | 0.92 | 0.86 | 0.60 |
  | top | 0.91 | 0.86 | 0.81 |

  The 3D reconstruction of the simulated holograms is consistent with that of the real ones,
  including the ring that appears when light enters through the flat face of the phantom. This
  run used the v1.1 index model (immersion correction); v1.2 uses the object's own index.
- Peak extraction matches the MATLAB extractor of the pipeline to 4·10⁻¹³ bins.
- **Reflection.** Backscatter probes object frequencies q_z = k(cos θi + cos θs), between 14.9 and
  29.8 rad/µm. At 12 px/µm, q·d reaches 2.5 rad per cell, close to Nyquist: the Yee grid resolves
  propagation well, but not the position of edges at that scale. Measured with Δn = 0.046, like
  the phantom:

  | Case (normal incidence) | 12 px/µm | 20 px/µm | 30 px/µm |
  |---|---|---|---|
  | Sharp sphere R = 1.5 µm vs Mie: complex error | 1.96 | 0.89 | 0.43 |
  | Sharp sphere: phase | ≈ π | −0.22 rad | −0.06 rad |
  | Sharp planar interface vs Fresnel: \|FDTD/Fresnel\| depending on the subcell position of the edge | 0.42–1.92 | 0.76–1.20 | — |
  | Planar interface: phase error | up to ±0.7 rad | up to ±0.4 rad | — |

  The solver is consistent: at 12 px/µm the FDTD result matches Born computed on the discrete Δε
  that Meep sees, and it converges to Mie when the grid is refined. The deviation comes from
  sampling the edges on the grid. The sphere is an extreme case, because its backscatter is the
  almost cancelled sum of two caps. In the same test, transmission has a 7–12 % error.

## Known limitations

- **`sf_calibration` factor and relative field (under review).** The empirical factor against Mie
  (1.19 at 12 px/µm) is applied only to E_s. But E_s and E_inc share Meep's normalization, so the
  ratio E_s/E_inc is already the physical one. Checked on the phantom slices: without the factor,
  the total intensity inside the object is 0.96–1.02 (as expected for a phase object); with it,
  1.66–1.77, and \|E_s/E_inc\| reaches 2.6 (it must be ≤ 2). The lateral slices no longer apply
  it. In `run` and `reflect`, S is inflated by 19 % at 12 px/µm.
- **Reflection at 12 px/µm is qualitative.** Where the interfaces reflect, the shape of the rings,
  the spectral support and the specular position are reliable. The amplitude of each interface is
  within a factor of ~0.4–2, and the phase within ±0.7 rad. In addition, the GT has a 0.121 µm
  voxel (Nyquist 25.9 rad/µm, below 2k = 29.8 rad/µm) and 1–2-voxel edges, which weaken the
  reflection near normal incidence compared with a real sharp interface. Reflection of the full
  phantom at 20 px/µm would need about 37 GB.
- Scalar field at the output (main component); the real polarization of the setup is unknown.
- The simulated holograms have no noise or speckle. The reference is uniform, and the
  illumination is an ideal plane wave, without the curvature of the real beam.
- Field scattered at large angles that leaves through the sides of the cell is lost. It depends
  on `lateral_margin_um` and `walkoff_margin_um`; see [RESOURCES.md](RESOURCES.md).
- With the v1.1 index model (n_immersion_gt), the central phase of the phantom was reproduced
  within 7 %. With v1.2, the accuracy depends on the index stored in the object.

## Modules

| Module | Content |
|---|---|
| `config.py` | Reading and validation of the run JSON |
| `acquisition.py` | `config_manip.txt`, PGM reading and Fourier geometry of the patch |
| `specular.py` | Specular peaks, FPS or list selection, and `real_subset/` |
| `illuminations.py` | Illumination lists: reading, generation and matching to real holograms |
| `objects.py` | Registered GT (medium index, automatic crop), h5 voxel format, view rotation, sphere |
| `scattered.py` | Discrete plane wave, cell geometry and scattered-field simulation |
| `optics.py` | NA pupil and angular-spectrum propagation |
| `camera.py` | Object → camera transform, demodulation, background, correlation, off-axis reference |
| `pipeline.py` | Commands `prepare`, `select`, `estimate`, `calibrate`, `run`, `compare` |
| `reflection.py` | Command `reflect`: backscatter, transfer to the opposite camera, reflection holograms |
| `free.py` | Command `free`: object + illumination list, no real data |
| `reference_views.py` | Command `reference`: lateral slices of the field in the FDTD domain and reference figures |
| `mie.py` | Exact angular spectrum (vector Mie) for validation |
| `example_sphere.py` | Self-contained example against Mie |
| `fdtd_meep_gui/` | PyQt5 interface (see [GUIDE.md](GUIDE.md)) |
