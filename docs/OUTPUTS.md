# Output files

**v1.2.0 · 2026-10-30 · Carlos Chacón**

What each file produced by the simulations contains (4Pi transmission and reflection, and free
mode), how to open it, and what it is for.

**Formats.** Simulated fields are saved as **`.npz`** (NumPy), not `.h5`. The only `.h5` file is
the **object** (the simulated index n(x,y,z)). The `.log` files are Meep's text output for each
illumination: they help with diagnosis but contain no data.

1. [Folder map](#1-folder-map)
2. [Formats and how to open them](#2-formats-and-how-to-open-them)
3. [Transmission: bottom/ and top/](#3-transmission-bottom-and-top)
4. [Reflection: refl_top/ and refl_bottom/](#4-reflection-refl_top-and-refl_bottom)
5. [Fields at the sensor (fields/*.npz)](#5-fields-at-the-sensor-fieldsnpz)
6. [Field at the reflection camera (camera_fields/*.npz)](#6-field-at-the-reflection-camera-camera_fieldsnpz)
7. [Object (.h5), calibration, reports and logs](#7-object-h5-calibration-reports-and-logs)
8. [Free mode: transmission/ and reflection/](#8-free-mode-transmission-and-reflection)
9. [What to use for what](#9-what-to-use-for-what)

---

## 1. Folder map

```
outputs/<run>/
├── object/object.h5             simulated object: index n(z,y,x)
├── calibration/                 normal-incidence simulations used for calibration (center_*.npz + .log)
│   ├── calibration.json         object → camera transform of each view
│   └── calibration.png          real vs simulated phase and amplitude of the central hologram
├── bottom/                      TRANSMISSION bottom (source on top, camera below)
│   ├── selection/               illumination_positions*.{csv,mat}, seleccion.txt, seleccion_k.png
│   ├── real_subset/             selected real holograms (i000…) + config
│   ├── sim_subset/              simulated holograms, same names + config + reports
│   ├── fields/                  Meep fields at the sensor (iNNN.npz) + Meep log (iNNN.log)
│   └── reference_views/         lateral slices of the field (optional, `reference`)
├── top/                         TRANSMISSION top (source below, camera on top), same as bottom/
├── refl_top/                    REFLECTION recorded by the top camera (bottom illuminations)
│   ├── fields/                  raw Meep reflection fields (iNNN.npz + iNNN.log)
│   ├── camera_fields/           exact complex spectrum of the +1 order (iNNN.npz)
│   ├── sim_subset/              reflection holograms without a parasitic specular
│   ├── sim_subset_spec0.01/     the same with a parasitic specular of amplitude 0.01
│   ├── refl_top_resumen.png     hologram, spectrum, |U| and phase for θ ≈ 0, 12, 25, 40°
│   └── refl_top_niveles.png     signal level and specular position vs θ
├── refl_bottom/                 REFLECTION recorded by the bottom camera (top illuminations), same
├── comparison/                  comparison.{json,png}: real vs simulated fields
└── reference_views/             fdtd_{zmin,zmax}_incidencia_normal.npz and 2 × 2 figures
```

**View mapping in the 4Pi**

| Folder | Source | Recording camera | Illumination directions |
|---|---|---|---|
| `bottom/` (transmission) | top | bottom | N selected bottom holograms |
| `top/` (transmission) | bottom | top | N selected top holograms |
| `refl_top/` | top | **top** | those of `bottom/` |
| `refl_bottom/` | bottom | **bottom** | those of `top/` |

The names `i000…` are the same in all folders of a view. For example,
`refl_top/sim_subset/i005.pgm` uses the same illumination as `bottom/sim_subset/i005.pgm`.
`seleccion.txt` maps each name to the original hologram of the acquisition.

## 2. Formats and how to open them

| Extension | Content | Python | MATLAB |
|---|---|---|---|
| `.pgm` | 8-bit hologram, 1024 × 1024 | `np.asarray(PIL.Image.open(f))` | `imread(f)` |
| `.npz` | Several named NumPy arrays (compressed) | `d = np.load(f); d['Es_Ex']` | Not directly: convert to `.mat` (see below) |
| `.h5` | Object volume + coordinates | `h5py.File(f)['n'][()]` → `n[z,y,x]` | `h5read(f,'/n')` → `n(x,y,z)` |
| `.csv` | Per-illumination tables | `csv` / `pandas.read_csv` | `readtable(f)` |
| `.json` | Parameters and calibrations | `json.load(open(f))` | `jsondecode(fileread(f))` |
| `.txt` | `config_manip.txt`, `recon.txt` (pipeline format), `seleccion.txt` | text | text |
| `.log` | Output of Meep or of a step | text | text |
| `.png` | Control figures | — | — |

To convert an `.npz` to `.mat`:

```python
import numpy as np
from scipy.io import savemat
d = np.load('refl_top/camera_fields/i000.npz')
savemat('i000.mat', {k: d[k] for k in d.files})
```

## 3. Transmission: bottom/ and top/

| File | What it is | Use |
|---|---|---|
| `real_subset/iNNN.pgm` | Selected **real** holograms, renamed i000… | Reference reconstruction with the same angles |
| `real_subset/Intensite_ref.pgm` | Real reference intensity | Pipeline |
| `real_subset/config_manip.txt` | Optical parameters of the camera (NA, N0, λ, focal lengths, CIRCLE_CX/CY…) with `NB_HOLO` = number selected | Pipeline; also used by the simulator |
| `real_subset/recon.txt` | Reconstruction parameters | Pipeline |
| `sim_subset/iNNN.pgm` | **Simulated** holograms, with the same names, orientation and carrier as the real ones | Go into the pipeline as they are |
| `sim_subset/Intensite_ref.pgm` | Synthetic, uniform reference | Pipeline |
| `sim_subset/sim_report.csv` | Per illumination: k in RAW and object conventions, simulated vs real specular peak, Δp (bins), saturation, Meep and wall time | Quality control |
| `sim_subset/sim_parameters.json` | Applied calibration (M, conj, t, focus), reference ratio, carrier, simulation parameters | Traceability |
| `fields/iNNN.npz` | Meep field at the sensor plane (see §5) | Direct Born/Rytov, other optics, reprocessing holograms |
| `fields/iNNN.log` | Meep output for that illumination | Diagnosis |
| `selection/illumination_positions_all.csv` | Specular peak of every real hologram: subpixel position, RAW k, θ, φ, quality, FPS rank | See the full scan |
| `selection/illumination_positions.{csv,mat}` | The selected ones + `new_filename` | Used by `run` and `reflect` |
| `selection/seleccion.txt` | new name ↔ original, kx, ky, θ | Traceability |
| `selection/seleccion_k.png` | (kx, ky) map of all and selected illuminations | Visual |

Columns of `illumination_positions*.csv`: `drow_subpix`/`dcol_subpix` are the specular position in
bins from the patch centre. `kx_norm_raw = dcol/REwald` and `ky_norm_raw = drow/REwald` are
normalized by k = 2π·n/λ. `theta_deg` is the angle in the medium.

## 4. Reflection: refl_top/ and refl_bottom/

| File | What it is | Use |
|---|---|---|
| `sim_subset/iNNN.pgm` | 8-bit reflection hologram **without a parasitic specular**: backscatter + reference only | What "should be seen" without parasitic reflections; a hard case for Rytov |
| `sim_subset_spec0.01/iNNN.pgm` | The same, plus a plane wave in the specular direction with amplitude 0.01 of the incident one | Mimics the parasitic reflection that creates a specular |
| `*/Intensite_ref.pgm` | Synthetic reference. \|R\| and the gray scale are one per view | Pipeline |
| `*/config_manip.txt`, `*/recon.txt` | Those of the **recording camera** (top for refl_top) | The pipeline uses the right geometry |
| `*/seleccion.txt` | Selection of the illumination view | Traceability |
| `*/sim_report.csv` | One row per hologram (columns below) | Control and analysis |
| `*/sim_parameters.json` | Illumination and camera views, both calibrations, a_ref, gray scale, parasitic r, simulation parameters | Traceability |
| `camera_fields/iNNN.npz` | Exact complex spectrum of the +1 order (see §6) | Noise-free reference for Rytov |
| `refl_*_resumen.png` | 4 angles: hologram, log\|spectrum\| with the expected specular (+) and the maximum (×), \|U\|, phase | Know what to look for in real data |
| `refl_*_niveles.png` | RMS and maximum of \|U\| vs θ; distance between the maximum and the specular, with and without parasitic | Summary of the view |

Columns of `sim_report.csv` (reflection):

| Column | Meaning |
|---|---|
| `new_filename`, `original`, `illumination_view` | Name, original hologram of the illumination, and the view it comes from |
| `theta_inc_deg` | Illumination angle in the medium |
| `kx/ky_norm_raw_illum` | k of the illumination, measured on the transmission camera |
| `kx/ky_norm_raw_specular`, `drow/dcol_specular` | Where the specular **should** be in the patch of the reflection camera (normalized k and bins) |
| `drow/dcol_peak` | Where the maximum of the patch of the 8-bit hologram is |
| `peak_to_specular_bins` | Distance between both: < 1 means the maximum marks the specular |
| `peak_snr` | Peak vs its surroundings |
| `parasitic_specular` | Amplitude r of the added specular (0 or 0.01) |
| `back_rms`, `back_max` | RMS and maximum of \|U\| at the camera, relative to \|E_inc\| (full 1024² FOV) |
| `saturated_fraction` | Fraction of pixels clipped at 255 |
| `meep_time`, `n_cells_M`, `wall_s` | Simulated time, cells (millions), wall time (s) |

In the raw reflection fields (`refl_*/fields/`) each `iNNN.npz` has the format of §5, with
`reflection = True` and the sensor below the object (z < 0).

## 5. Fields at the sensor (fields/*.npz)

This is what Meep computes, before any optics. There is one file per illumination, in the frame
of the simulated object. For the phantom it takes about 20 MB.

| Variable | Type / shape | Meaning |
|---|---|---|
| `Es_Ex`, `Es_Ey`, `Es_Ez` | complex [ny, nx] | **Field scattered** by the object at the sensor plane, phasor at the light frequency (exp(−iωt) convention) |
| `Einc_Ex`, `Einc_Ey`, `Einc_Ez` | complex [ny, nx] | Discrete incident wave, evaluated on the same plane with the same normalization |
| `x`, `y` | µm, 1D | Plane coordinates (step 1/resolution µm) |
| `z` | µm | Height of the plane: > 0 in transmission (above the object), < 0 in reflection (below) |
| `kx`, `ky`, `kappa` | rad/µm | Incident wave vector. `kappa` comes from the discrete Yee dispersion |
| `kx_norm`, `ky_norm` | — | Transverse k / (2π·n/λ) in the object frame |
| `pol` | complex (3,) | Incident polarization (≈ x, projected ⟂ k) |
| `orientation` | text | Face through which light enters: `zmin` (object as is) or `zmax` (rotated 180° about y). Older files say `bottom`/`top`, equivalent to `zmin`/`zmax` |
| `reflection` | bool | Reflection only |
| `cell`, `cell_center`, `walkoff`, `n_cells` | µm / cells | Meep cell geometry |
| `meep_time`, `resolution` | — | Simulated time and px/µm |

The useful quantity is the **relative field** S = Es/Einc of the main component (`Ex`):

- **Transmission:** the total field is E_inc·(1 + S). Rytov gives ln(1 + S) and Born, S.
- **Reflection:** there is no incident wave on that plane. S is the backscatter with the ramp of
  the incident transverse k already removed (specular reflection preserves that k).

```python
d = np.load('bottom/fields/i010.npz')
S = d['Es_Ex'] / d['Einc_Ex']            # [y, x], in the object frame
```

These fields are at the sensor plane, without NA pupil or focus, and in object coordinates. To
carry them to the camera, use the calibration and `fdtd_meep.optics.focus_field` +
`fdtd_meep.camera.to_camera_grid`, as in `pipeline.run_view` and `reflection.camera_field`.

### Lateral slices (reference_views/fdtd_*_incidencia_normal.npz)

The variables are those of the table above, plus the slices of the central illumination:
`slice_planes`; `slice_<XZ|YZ>_{x,y,z}` (µm, simulation frame); `slice_<plane>_Es_{Ex,Ey,Ez}` and
`slice_<plane>_Einc_{Ex,Ey,Ez}`, of shape [nx, ny, nz] with one axis of size 1, and `z_pad_um`.
Here E_s does **not** carry the `sf_calibration` factor, so E_s + E_inc is directly the total
field.

## 6. Field at the reflection camera (camera_fields/*.npz)

This is the cross term U·R* of the reflection hologram, already on the camera, without
quantization or the autocorrelation term. It is cropped to the +1-order patch that the pipeline
uses. The numbers in parentheses are those of the phantom acquisition.

| Variable | Meaning |
|---|---|
| `spectrum_patch` | complex [304, 304]: +1-order spectrum of the patch, normalized by \|R\|, without a parasitic specular |
| `patch_center` | Patch centre, 1-based (153), where k_camera = 0 |
| `drow_specular`, `dcol_specular` | Expected specular position in bins from the centre: 0-based index = `patch_center − 1 + d` |
| `kx/ky_norm_raw_specular` | The same, in normalized k (RAW convention) |
| `fmaxHolo`, `REwald` | Radius of the NA pupil in bins (152) and Ewald radius (174.8 bins) |
| `dx_patch_um` | Pixel of the field reconstructed at patch resolution (0.2425 µm; FOV 73.7 µm) |
| `a_ref`, `gray_scale` | Reference amplitude and gray scale used for the `.pgm` files |

Complex camera field at patch resolution, relative to the incident amplitude. It carries the
specular carrier, that is, the tilt of the illumination:

```python
d = np.load('refl_top/camera_fields/i000.npz')
P = d['spectrum_patch']; n = P.shape[0]
U = np.fft.ifft2(np.fft.ifftshift(P)) * n * n        # |U| ≈ amplitude / |E_inc|
c = int(d['patch_center']) - 1
spec = (c + d['drow_specular'], c + d['dcol_specular'])   # where the specular would be in P
```

For a parasitic specular of amplitude r, add r·exp(i k·r) to U, that is, a peak at `spec`: other
amplitudes can be tested without simulating again.

## 7. Object (.h5), calibration, reports and logs

**Object (`object/object.h5`).** The index that is simulated.

| Content | Detail |
|---|---|
| `/n` | float32; h5py reads it as `n[z, y, x]` (phantom: 116 × 266 × 210), MATLAB as `n(x,y,z)` |
| `/coords/{x,y,z}_um` | Voxel centres in µm, centred at 0 (phantom: 0.1213 µm voxel, 25.5 × 32.3 × 14.1 µm volume) |
| attribute `metadata_json` | Origin (registered GT), crop used, voxel, n_medium, Δn at the crop border, n min/max |

In v1.2 the index is n = n_medium + Δn_GT, with n_medium taken from the object and no immersion
correction (see [DOCUMENTATION.md](DOCUMENTATION.md#1-object-prepare)).

**Calibration.**

| File | Content |
|---|---|
| `calibration/calibration.json` | Per view: the best match (`best`) between the simulation and the central real hologram. It includes the symmetry M, `conj` (−1 order), the shift `t_um`, the focal plane `z_focus_um`, the orientation and the correlations; `all` keeps every candidate tried |
| `calibration/center_{zmin,zmax}.npz` | Normal-incidence simulations, in the format of §5 |
| `calibration/calibration.png` | Real vs simulated phase and amplitude of the central hologram |

**Comparison.** `comparison/comparison.json` has the complex correlation between the real and
simulated demodulated fields (view, hologram, θ). `comparison.png` shows it.

**Logs.**

| File | Content |
|---|---|
| `iNNN.log` (next to each `.npz`) | Meep output: MPI processes, cell size, `set_epsilon` time, one step every ~1 s and the field decay; at the end, `run 0 finished at t = … (N timesteps)` and the total time. Useful if an illumination fails, takes too long (swapping: s/step ≫ 0.15) or does not converge |
| GUI console | One line per illumination: name, original, θ, cells, memory and, in transmission, Δp |

## 8. Free mode: transmission/ and reflection/

```
outputs/<run>/
├── object/object.h5
├── transmission/
│   ├── fields/<name>.npz         field at the sensor (format of §5)
│   ├── free_report.csv           per illumination: kx, ky, θ, φ, RMS and max |S|, sensor z, cells, times
│   ├── free_parameters.json      mode, list, component, n_medium, λ, object and simulation parameters
│   ├── illuminations_used.csv    the list actually simulated (kx, ky, θ, φ)
│   └── transmission_resumen.png  |1 + S| and phase at the sensor for up to 4 illuminations
├── reflection/                   the same; |S| (backscatter only) and phase
└── mie_check/                    optional: solver check against Mie (sphere_vs_mie.png, report.json)
```

## 9. What to use for what

| Goal | What to use |
|---|---|
| Reconstruct with the pipeline as if it were real data | `*/sim_subset/` (or `sim_subset_spec0.01/`): holograms + `config_manip.txt` + `recon.txt` + `Intensite_ref.pgm` |
| Compare with real data at the same angles | `bottom|top/real_subset/` vs `bottom|top/sim_subset/` |
| Know what to look for in real reflection holograms | `refl_*_resumen.png`, `refl_*_niveles.png` |
| Develop a Rytov normalization without a specular | `refl_top/camera_fields/` (noise-free truth) vs `refl_top/sim_subset/` (8 bits). `drow/dcol_specular` tell where the specular should be |
| Test the effect of a parasitic specular | `sim_subset_spec0.01/`, or add r to `spectrum_patch` (§6) |
| Test Born/Rytov without camera optics | `fields/`: S = Es_Ex/Einc_Ex |
| Reference figure of the propagation (paper, talk) | `reference_views/lateral_reference_4casos_*.{png,svg}` or `<view>/reference_views/` |
| Simulate another object | The object `.h5` format (§7); see `matlab/exportar_objeto_meep.m` |

**Reflection limitation.** At 12 px/µm reflection is **qualitative**. The position of the
interfaces, the shape of the rings, the spectral support and the specular position are reliable.
The amplitude of each interface can vary by a factor of ~0.4 to 2, and the phase by ±0.7 rad.
See the validation in [DOCUMENTATION.md](DOCUMENTATION.md#validation).
