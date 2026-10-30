# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón
# Free mode: any object + an illumination list, no real holograms or calibration needed.

import csv
import json
import time

import numpy as np

from .illuminations import read_illumination_list, write_illumination_list

MODES = ('transmission', 'reflection')


def run_free(pipe, mode, limit=None):
    """Simula cada iluminación de la lista con el objeto tal cual (orientación zmin).

    Las direcciones (kx, ky) están en el marco del objeto. Salidas en outputs/<run>/<modo>/:
    fields/<name>.npz (mismo formato que run/reflect), free_report.csv,
    illuminations_used.csv y <modo>_resumen.png.
    """
    if mode not in MODES:
        raise ValueError(f'Modo desconocido: {mode} (transmission o reflection).')
    reflection = mode == 'reflection'
    rows = read_illumination_list(pipe.cfg.illuminations_file())
    if any('kx' not in r for r in rows):
        raise ValueError('En modo libre cada iluminación necesita kx,ky o theta_deg,phi_deg.')
    rows = rows[:limit] if limit else rows
    vdir = pipe.view_dir(mode)
    fields_dir = vdir / 'fields'
    fields_dir.mkdir(parents=True, exist_ok=True)
    write_illumination_list(vdir / 'illuminations_used.csv', rows)
    comp = pipe.cfg.simulation.polarization
    records = []
    for j, row in enumerate(rows):
        npz = fields_dir / f"{row['name']}.npz"
        t0 = time.time()
        if not npz.exists():
            lay = pipe.layout('zmin', row['kx'], row['ky'], reflection)
            pipe.check_memory(lay)
            print(f"[{mode} {j + 1}/{len(rows)}] {row['name']} θ={row['theta_deg']:.1f}° "
                  f"φ={row['phi_deg']:.1f}° celdas={lay['n_cells'] / 1e6:.1f} M "
                  f"(~{lay['memory_gb']:.1f} GB)", flush=True)
            pipe.run_sim('zmin', row['kx'], row['ky'], npz, reflection=reflection)
        S, x, y, d = pipe.relative_field(npz, comp)
        records.append(dict(name=row['name'], kx=row['kx'], ky=row['ky'],
                            theta_deg=row['theta_deg'], phi_deg=row['phi_deg'],
                            rms_S=float(np.sqrt(np.mean(np.abs(S) ** 2))),
                            max_S=float(np.abs(S).max()), sensor_z_um=float(d['z']),
                            n_cells_M=float(d['n_cells']) / 1e6, meep_time=float(d['meep_time']),
                            wall_s=time.time() - t0))
    with open(vdir / 'free_report.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader()
        w.writerows(records)
    (vdir / 'free_parameters.json').write_text(json.dumps(dict(
        mode=mode, illuminations=str(pipe.cfg.illuminations_file()), component=comp,
        n_medium=pipe.n_medium, wavelength_um=pipe.wavelength, object=vars(pipe.cfg.object),
        simulation=vars(pipe.cfg.simulation)), indent=2, default=str))
    figure(vdir / f'{mode}_resumen.png', pipe, mode, rows, fields_dir, comp)
    return records


def figure(path, pipe, mode, rows, fields_dir, comp, max_cols=4):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    picks = np.unique(np.linspace(0, len(rows) - 1, min(max_cols, len(rows))).round().astype(int))
    fig, ax = plt.subplots(2, len(picks), figsize=(4.2 * len(picks), 8.2), squeeze=False)
    for col, j in enumerate(picks):
        row = rows[j]
        S, x, y, d = pipe.relative_field(fields_dir / f"{row['name']}.npz", comp)
        # Transmisión: campo total relativo 1 + S; reflexión: solo la retrodispersión S
        U = S if mode == 'reflection' else 1 + S
        e = [x[0], x[-1], y[0], y[-1]]
        lab = '|S|' if mode == 'reflection' else '|1 + S|'
        im = ax[0, col].imshow(np.abs(U), origin='lower', extent=e, cmap='hot')
        fig.colorbar(im, ax=ax[0, col], fraction=0.046)
        ax[0, col].set_title(f"{mode} {row['name']} θ={row['theta_deg']:.0f}° φ={row['phi_deg']:.0f}°\n"
                             f"{lab} en el sensor (z = {float(d['z']):+.2f} µm)", fontsize=8)
        im = ax[1, col].imshow(np.angle(U), origin='lower', extent=e, cmap='twilight',
                               vmin=-np.pi, vmax=np.pi)
        fig.colorbar(im, ax=ax[1, col], fraction=0.046)
        ax[1, col].set_title('fase (rad)', fontsize=8)
        for a in ax[:, col]:
            a.set_xlabel('x (µm)', fontsize=7)
            a.set_ylabel('y (µm)', fontsize=7)
            a.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)
