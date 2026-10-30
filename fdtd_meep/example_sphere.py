# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón

import json
import math
import time
from pathlib import Path

import numpy as np

from .config import SimulationParams
from .mie import angular_spectrum
from .objects import sphere_volume, write_voxel_h5
from .scattered import simulate_scattered


def yee_kz(KX, KY, pw, n_medium):
    d, dt, w = pw['d'], pw['dt'], pw['omega']
    s2 = (d / dt) ** 2 * n_medium ** 2 * math.sin(w * dt / 2) ** 2 \
        - np.sin(KX * d / 2) ** 2 - np.sin(KY * d / 2) ** 2
    ok = (s2 > 0) & (s2 < 1)
    return np.where(ok, 2 / d * np.arcsin(np.sqrt(np.clip(s2, 0, 1))), np.nan)


def run_case(n, step, n_sphere, n_medium, radius, wavelength, kx_norm, ky_norm, params, sensor_z,
             reflection=False):
    t0 = time.time()
    r = simulate_scattered(n, step, (0.0, 0.0, 0.0), n_medium, wavelength, kx_norm, ky_norm,
                           sensor_z, params)
    pw, x, y, z = r['pw'], r['x'], r['y'], r['z']
    E = r['fields']['Ex'] / r['S']
    d = x[1] - x[0]
    pad = 4
    F = np.fft.fft2(E, s=(len(y) * pad, len(x) * pad)) * d * d
    KX, KY = np.meshgrid(2 * np.pi * np.fft.fftfreq(len(x) * pad, d),
                         2 * np.pi * np.fft.fftfreq(len(y) * pad, d))
    F *= np.exp(-1j * (KX * x[0] + KY * y[0]))
    sz = 1 if reflection else -1
    A_fd = F * np.exp(sz * 1j * yee_kz(KX, KY, pw, n_medium) * z) / np.cos(KY * d / 2)
    k = 2 * np.pi * n_medium / wavelength
    k_inc = [pw['kx'], pw['ky'], math.sqrt(k ** 2 - pw['kx'] ** 2 - pw['ky'] ** 2)]
    A_mie = angular_spectrum(KX, KY, k_inc, pw['pol'], n_sphere, n_medium, radius, wavelength,
                             backward=reflection)[..., 0]
    region = (np.hypot(KX, KY) < 0.8 * k) & np.isfinite(A_fd)
    m = region & (np.abs(A_mie) > 0.05 * np.abs(A_mie[region]).max())
    a, b = A_fd[m], A_mie[m]
    ratio = complex(np.vdot(a, b) / np.vdot(a, a))
    return dict(kx_norm=kx_norm, ky_norm=ky_norm,
                theta_deg=math.degrees(math.asin(math.hypot(kx_norm, ky_norm))),
                rel_error=float(np.linalg.norm(b - a) / np.linalg.norm(b)),
                amplitude_ratio_mie_over_fdtd=abs(ratio), phase_offset_rad=float(np.angle(ratio)),
                meep_time=r['meep_time'], wall_s=time.time() - t0,
                cells_M=r['layout']['n_cells'] / 1e6), (np.fft.fftshift(A_fd), np.fft.fftshift(A_mie),
                                                        np.fft.fftshift(m))


def run(output, resolution=12, angles=((0.0, 0.0), (0.3, 0.0), (0.45, 0.35)),
        radius=1.5, dn=0.03, n_medium=1.5, wavelength=0.6328, reflection=False, voxel_um=0.1,
        edge_um=0.05, pml_far_um=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    vox = voxel_um
    n = sphere_volume(radius, dn, n_medium, vox, edge_um=edge_um)
    write_voxel_h5(output / 'sphere.h5', n, vox, dict(radius_um=radius, dn=dn, n_medium=n_medium,
                                                       edge_um=edge_um))
    params = SimulationParams(resolution=resolution, lateral_margin_um=3.0, mpi_processes=1,
                              pml_far_um=pml_far_um)
    sensor_z = (-1 if reflection else 1) * (n.shape[0] * vox / 2 + 0.4)
    report, maps = [], []
    for kx, ky in angles:
        res, m = run_case(n, [vox] * 3, n_medium + dn, n_medium, radius, wavelength, kx, ky,
                          params, sensor_z, reflection)
        print(f"θ={res['theta_deg']:5.1f}°  error relativo frente a Mie = {res['rel_error']:.3f}  "
              f"fase = {res['phase_offset_rad']:+.3f} rad  ({res['wall_s']:.0f} s)", flush=True)
        report.append(res)
        maps.append(m)
    (output / 'report.json').write_text(json.dumps(dict(
        resolution=resolution, radius_um=radius, dn=dn, n_medium=n_medium,
        wavelength_um=wavelength, reflection=reflection, voxel_um=vox, edge_um=edge_um,
        pml_far_um=pml_far_um, cases=report), indent=2))
    figure(output / 'sphere_vs_mie.png', report, maps)
    return report


def figure(path, report, maps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(len(maps), 3, figsize=(13, 4 * len(maps)), squeeze=False)
    for i, ((A_fd, A_mie, m), r) in enumerate(zip(maps, report)):
        vmax = np.abs(np.where(m, A_mie, 0)).max()
        pr, pc = np.unravel_index(int(np.argmax(np.abs(np.where(m, A_mie, 0)))), A_mie.shape)
        h = min(A_fd.shape) // 6
        c = (slice(max(pr - h, 0), pr + h), slice(max(pc - h, 0), pc + h))
        ax[i, 0].imshow(np.abs(np.where(m, A_fd, 0))[c], vmax=vmax, cmap='magma')
        ax[i, 0].set_title(f"|A| FDTD (θ={r['theta_deg']:.0f}°)")
        ax[i, 1].imshow(np.abs(np.where(m, A_mie, 0))[c], vmax=vmax, cmap='magma')
        ax[i, 1].set_title('|A| Mie exacto')
        ax[i, 2].imshow(np.where(m, np.angle(A_fd * np.conj(A_mie)), np.nan)[c],
                        vmin=-0.5, vmax=0.5, cmap='RdBu_r')
        ax[i, 2].set_title(f"Δfase (rad), error relativo {r['rel_error']:.3f}")
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)
