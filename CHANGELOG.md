# Changelog

## 1.2.0 — 2026-10-30

- **PyQt5 interface** (`fdtd_meep_gui`, `./run_gui.sh`): parameters, step queue with progress
  and stop, results viewer (object, illuminations, holograms, sensor field, report), saved
  figures and console.
- **Free mode** (`python -m fdtd_meep free`): any object plus an illumination list, in
  transmission and/or reflection, without real holograms or calibration.
- **Illumination lists** (CSV with `kx,ky`, `theta_deg,phi_deg` or `filename`), usable for the
  4Pi selection (`selection.lists`) and in free mode (`illuminations`); a list generator (spiral
  or rings).
- Sphere object source (`object.source = "sphere"`).
- **Medium index from the object, no index correction**: n = n_medium + Δn, with n_medium taken
  from the object file (`.mat` `meta.n_medium` / `meta.phantom_n_immersion`, `.h5` attribute
  `n_medium`) unless set in the config. The v1.1 immersion correction (`n_immersion_gt`) is
  applied only if requested explicitly. `prepare` warns when n_medium differs from N0.
- **Automatic crop** of the registered GT (support + `crop_pad_um`) when `crop_xyz` is empty,
  and a warning if the object touches the crop border.
- `install.sh`, an English README and `docs/GUIDE.md`.
- The solver, calibration, transmission and reflection are unchanged. Phantom results differ
  from 1.1.0 only through the object index (no immersion correction).

## 1.1.0 — 2026-08-26

- 4Pi reflection (`reflect`), lateral reference views (`reference`).

## 1.0.0 — 2025-07-05

- Scattered-field FDTD core, selection, calibration, transmission, comparison, Mie validation.
