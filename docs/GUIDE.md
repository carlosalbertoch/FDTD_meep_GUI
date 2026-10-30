# FDTD_meep GUI — User guide

**v1.2.0 · 2026-10-30 · Carlos Chacón**

The interface writes a JSON configuration and runs the same commands as the command line
(`python -m fdtd_meep <step> --config <json>`), one after the other, in a separate process.
Closing the window does not change any result already written.

## Left panel: parameters

### Project

| Field | What it does |
|---|---|
| Workflow | *4Pi acquisition*: real holograms of both views. *Free simulation*: any object plus an illumination list, with no real data |
| Config | The JSON file being edited. **Open…** loads one, **Save** writes the form to it, **Save as…** writes a copy. **Run** always saves first |
| Run name | Name of the run (`name` in the JSON) |
| Output folder | Where the results go (default `outputs/<run name>`). Point it at an earlier run to view its results |

### Input holograms (4Pi only)

One folder per view, each with `config_manip.txt`, `Intensite_ref.pgm` and `iNNN.pgm`. Under the
folders the GUI shows λ, N0, NA, the image size and the number of holograms read from
`config_manip.txt`. If the folder is not valid, it shows a warning instead.

### Object

| Source | Fields | Notes |
|---|---|---|
| `sphere` | Radius, Δn, voxel | Smooth edge (tanh, 0.05 µm). Best for tests and validation |
| `voxel_h5` | File | Your own volume. Export it from MATLAB with `matlab/exportar_objeto_meep.m` (format in [DOCUMENTATION.md](DOCUMENTATION.md#your-own-objects-from-matlab)) |
| `registered_gt` | File (.mat v7.3), crop, crop margin | Registered ground truth of the phantom (Δn), cropped to a box around the object |

**Medium index.** No index correction is applied: n = n_medium + Δn. With *n medium* on *from
object file*, the index stored in the object is used (`.mat`: `meta.n_medium`, otherwise
`meta.phantom_n_immersion`; `.h5`: attribute `n_medium`). Type a value to override it. The line
under the object shows the value in use, where it comes from, and a warning if it differs from
N0 of `config_manip.txt`. N0 still sets the camera optics. A sphere has no file, so it takes
*n medium*, or N0 in 4Pi mode.

**Crop.** The GT is a large cube (608³ voxels, 74 µm); only a box around the object is
simulated. The box is given by its faces along each axis, in GT voxels: `[[x0,x1],[y0,y1],[z0,z1]]`,
0-based, with the end excluded (Python slicing). In MATLAB, with V(y, x, z), this is
`V(y0+1:y1, x0+1:x1, z0+1:z1)`. Leave the field empty, or click **Auto**, to use the box of the
support (|Δn| > threshold) plus *Crop margin* of medium on every side. The cropped box is
centred at (0, 0, 0); *Calibrate* finds the shift, symmetry and focus with respect to each
camera. A warning appears if the object touches the border of the box.

**4Pi flip.** It is automatic: the object is simulated as is (*zmin*, light enters through the
z-min face) and rotated 180° about y (*zmax*). *Calibrate* assigns each view to one face, and the
two views always enter through opposite faces. *Reflect* reuses that assignment.

*Wavelength* is shown only in free mode; in 4Pi mode λ comes from the acquisition.

### Simulation

| Field | Default | Meaning |
|---|---|---|
| Resolution | 12 px/µm | Calibrated at 12, 14 and 20. Below 12 there are fewer than 5 points per wavelength in n = 1.55 |
| Polarization | Ex | Ex or Ey in the object frame, projected perpendicular to k |
| PML | 1.0 µm | Absorbing layer thickness |
| Lateral / walk-off / z margin | 0.8 / 1.0 / 0.4 µm | Space around the object. Oblique light leaves the object shifted sideways (walk-off), and the cell grows only toward that side |
| Sensor gap | 0.4 µm | Distance from the object box to the sensor plane |
| MPI processes | 4 | Processes per illumination. On an M-series Mac with 8 performance cores, 4 is optimal |
| Max. memory | no limit | An illumination whose estimated memory exceeds this limit is not launched |

### Illuminations

- **4Pi, Automatic (FPS)**: farthest-point sampling in (kx, ky) of the real specular peaks. It
  starts at the most central direction and takes N per view.
- **4Pi, From list file**: one CSV list per view. Each direction is matched to the nearest real
  hologram; with a `filename` column, the named hologram is used.
- **Free**: one CSV list. **Generate list…** creates N directions inside a cone θ ≤ θmax, as a
  Fibonacci spiral (uniform coverage) or concentric rings.

The format is described in the [README](../README.md#illumination-list-format). Below the list,
the GUI shows the number of illuminations and the maximum angle.

### Process

1. Tick *Transmission* and/or *Reflection*. The step list is rebuilt from that choice.
2. Untick the steps that are already done. Steps can also be resumed: fields that already exist
   are not recomputed.
3. *Illuminations per view* limits the run to the first N illuminations, which is useful for a
   quick test (`--limit`).
4. **▶ Run selected steps** saves the config and starts the run. **■ Stop** ends the whole
   process group, including mpirun and its ranks. An illumination that was cut off leaves no
   field, so it is simply recomputed on the next run.
5. **Open output folder** opens the results in Finder (or your file manager).

## Right panel: results

### Results viewer

The top bar selects the **View** (bottom, top, refl_top, refl_bottom; or transmission and
reflection in free mode), the **Illumination**, the **Hologram** source (`input` = acquisition
folder, `real` = selected real subset, `sim` = simulated, `sim_spec<r>` = simulated with a
parasitic specular), the **Colormap**, and **Log scale**.

| Panel | Content |
|---|---|
| Material n | Central slice XZ / YZ / XY of the prepared object, or a preview of the input object before *Prepare* |
| Illuminations | Directions in (kx/k, ky/k): all detected peaks in grey, the selected ones in color. The current illumination is circled. The circle marks kx² + ky² = 1 |
| Hologram | *Intensity* (gray levels), and from the demodulated +1 order: *φ wrapped*, *φ unwrapped* and *Amplitude*. The phase background is removed with the same polynomial fit as `calibrate` |
| Sensor field amplitude / phase | Relative field at the sensor plane from `fields/<name>.npz`: \|1 + S\| in transmission, \|S\| (backscatter only) in reflection, with S = E_s/E_inc |
| Details | Report of the illumination (θ, k, Δp or distance to the specular, cells, Meep time, wall time) and the calibration of each view |

### Saved figures

Lists every PNG in the output folder: the selection coverage, the calibration, the comparison,
the transmission and reflection summaries, the lateral reference views, and the Mie check. Use
the toolbar to zoom. Double-click a figure to open it in the system viewer.

### Console

Shows the full output of every command, prefixed with the command line it runs. When a step
fails, its error is here. Copy that command to run the step by hand.

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `View bottom: … has no config_manip.txt` | The view folder must contain the acquisition files directly, not a parent folder |
| `La celda estimada … supera max_memory_gb` | Raise *Max. memory*, reduce the margins, or lower the resolution. Run *Estimate memory* first |
| `La iluminación no es propagante en la malla discreta` | kx² + ky² is too close to 1 for this resolution: reduce the angle or increase the resolution |
| A step takes 10× longer than expected | The machine is swapping. Close MATLAB and browsers, or set *MPI processes* to 2 |
| `calibrate`/`run` fails on a fresh run | The steps depend on each other: *Prepare → Select → Calibrate → Transmission/Reflection*. Run them in this order at least once |
| Nothing in *Sensor field* | That view or illumination has not been simulated yet, or *Output folder* points to another run |
