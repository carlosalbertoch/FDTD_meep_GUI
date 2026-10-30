# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón

"""Figuras laterales de referencia: el campo dentro del dominio FDTD (cortes XZ / YZ) para
la iluminación central (kx = ky = 0), en transmisión y reflexión.

Salen de monitores DFT de Meep sobre el corte, no de los campos del sensor ni de
camera_fields. Las figuras se dibujan en el marco del laboratorio (el del objeto tal cual,
orientación zmin): +z es la parte superior del microscopio en todas las vistas, así que el
phantom aparece igual y lo que cambia es el sentido de la luz y la cámara que registra.
Una vista de transmisión y la vista de reflexión con la misma iluminación comparten la
misma simulación: solo cambia la cámara.
"""

import dataclasses
import math

import numpy as np

from .objects import read_voxel_h5
from .scattered import cell_layout, discrete_plane_wave, volume_axes

TITLES = {'top': 'Transmisión top', 'bottom': 'Transmisión bottom',
          'refl_top': 'Reflexión top', 'refl_bottom': 'Reflexión bottom'}
COMPS = ('Ex', 'Ey', 'Ez')


def case_geometry(pipe, view):
    cal = pipe.calibration()
    if view in pipe.cfg.reflection.views:
        spec = pipe.cfg.reflection.views[view]
        orientation = cal[spec['illumination']]['best']['orientation']
        return dict(kind='reflection', orientation=orientation, camera=spec['camera'],
                    light_up=orientation == 'zmin')
    orientation = cal[view]['best']['orientation']
    return dict(kind='transmission', orientation=orientation, camera=view,
                light_up=orientation == 'zmin')


def slice_center_sim(pipe, orientation):
    rv = pipe.cfg.reference_views
    c = np.array(pipe.cfg.object.center_um, float)
    if rv.position_um is not None:
        c[0] = c[1] = rv.position_um
    if orientation == 'zmax':
        c[0], c[2] = -c[0], -c[2]
    return c


def reference_npz(pipe, orientation):
    return pipe.out / 'reference_views' / f'fdtd_{orientation}_incidencia_normal.npz'


def simulate_reference(pipe, orientation):
    rv = pipe.cfg.reference_views
    npz = reference_npz(pipe, orientation)
    if npz.exists():
        return npz
    p = dataclasses.replace(pipe.cfg.simulation, z_margin_um=pipe.cfg.simulation.z_margin_um + rv.z_pad_um)
    n, step, center, sensor_z = pipe.load_object(orientation)
    pw = discrete_plane_wave(p.resolution, pipe.wavelength, pipe.n_medium, 0.0, 0.0, p.polarization)
    lay = cell_layout(n.shape, step, center, sensor_z, pw, p)
    pipe.check_memory(lay)
    print(f"Referencia lateral, entrada por {orientation}: celdas={lay['n_cells'] / 1e6:.1f} M "
          f"(~{lay['memory_gb']:.1f} GB)", flush=True)
    c = slice_center_sim(pipe, orientation)
    pipe.run_sim(orientation, 0.0, 0.0, npz, extra_args=[
        '--slices', ','.join(rv.planes), f'--z-pad={rv.z_pad_um}', f'--slice-center={c[0]},{c[1]},{c[2]}'])
    return npz


def lab_slice(d, plane, orientation):
    """Corte en el marco del laboratorio: arrays [z, h] con h = x (XZ) o y (YZ)."""
    h_axis = 'x' if plane == 'XZ' else 'y'
    h, z = d[f'slice_{plane}_{h_axis}'], d[f'slice_{plane}_z']
    sq = lambda a: (a[:, 0, :] if plane == 'XZ' else a[0, :, :]).T
    Es = {c: sq(d[f'slice_{plane}_Es_{c}']) for c in COMPS}
    Ei = {c: sq(d[f'slice_{plane}_Einc_{c}']) for c in COMPS}
    if orientation == 'zmax':
        z = -z[::-1]
        Es = {c: v[::-1] for c, v in Es.items()}
        Ei = {c: v[::-1] for c, v in Ei.items()}
        if plane == 'XZ':
            h = -h[::-1]
            Es = {c: v[:, ::-1] for c, v in Es.items()}
            Ei = {c: v[:, ::-1] for c, v in Ei.items()}
    other = d[f'slice_{plane}_{"y" if plane == "XZ" else "x"}'][0]
    if orientation == 'zmax' and plane == 'YZ':
        other = -other
    return h, z, Es, Ei, float(other)


def phantom_slice(pipe, plane, h, z, pos):
    from scipy.interpolate import RegularGridInterpolator
    n, step = read_voxel_h5(pipe.object_path())
    ax, ay, az = volume_axes(n.shape, step, pipe.cfg.object.center_um)
    f = RegularGridInterpolator((az, ay, ax), n, bounds_error=False, fill_value=pipe.n_medium)
    Z, H = np.meshgrid(z, h, indexing='ij')
    P = np.full(Z.shape, pos)
    pts = np.stack([Z, P, H] if plane == 'XZ' else [Z, H, P], -1)
    return f(pts.reshape(-1, 3)).reshape(Z.shape)


def compute_display_intensity(I, mode, alpha):
    if mode == 'raw':
        return I / np.percentile(I, 99.9)
    vmax = np.percentile(I, 99.9)
    return np.clip(np.arcsinh(alpha * I) / np.arcsinh(alpha * vmax), 0, 1)


def auto_alpha(I, n2d, z, n_medium, camera_top):
    """α = 1/I típica en el lado de la cámara, fuera del objeto: lineal hasta ese nivel,
    logarítmica por encima. Así se ve lo que llega a la cámara sin saturar el resto."""
    inside = n2d > n_medium + 0.005
    rows = np.flatnonzero(inside.any(1))
    zone = (z > z[rows[-1]]) if camera_top else (z < z[rows[0]])
    m = zone[:, None] & ~inside
    ref = np.median(I[m]) if m.any() else np.median(I)
    return 1.0 / max(ref, 1e-12)


def overlay_phantom_contour(ax, h, z, n2d, n_medium, color='white'):
    level = n_medium + 0.5 * (n2d.max() - n_medium)
    ax.contour(h, z, n2d, levels=[level], colors=color, linewidths=0.6, alpha=0.85)


def plot_lateral_field(ax, D, h, z, n2d, n_medium, geom, rv, title, alpha):
    import matplotlib.patheffects as pe
    im = ax.imshow(D, origin='lower', extent=[h[0], h[-1], z[0], z[-1]], cmap=rv.colormap,
                   vmin=0, vmax=1, interpolation='bilinear', aspect='equal')
    im.set_rasterized(True)
    if rv.overlay_contour:
        overlay_phantom_contour(ax, h, z, n2d, n_medium)
    stroke = [pe.withStroke(linewidth=2, foreground='black')]
    y0, y1 = (0.06, 0.24) if geom['light_up'] else (0.94, 0.76)
    ax.annotate('', xy=(0.06, y1), xytext=(0.06, y0), xycoords='axes fraction',
                arrowprops=dict(arrowstyle='-|>', color='white', lw=1.4))
    ax.text(0.09, (y0 + y1) / 2, r'$k_{inc}$', transform=ax.transAxes, color='white', fontsize=9,
            va='center', path_effects=stroke)
    cam_top = geom['camera'] == 'top'
    ax.text(0.5, 0.975 if cam_top else 0.025, f"cámara {geom['camera']}", transform=ax.transAxes,
            color='white', fontsize=8, ha='center', va='top' if cam_top else 'bottom', path_effects=stroke)
    ax.set_xlabel(f"{'x' if 'XZ' in title else 'y'} (µm)", fontsize=8)
    ax.set_ylabel('z (µm)  ↑ parte superior del microscopio', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=9)
    return im


def colorbar_label(rv, alpha):
    if rv.intensity_mode == 'raw':
        return r'$I/I_{99.9\%}$,  $I=|E|^2/|E_{inc}|^2$'
    return rf'arcsinh($\alpha I$)/arcsinh($\alpha I_{{máx}}$),  $I=|E|^2/|E_{{inc}}|^2$,  $\alpha$={alpha:.3g}'


def case_images(pipe, view, plane):
    geom = case_geometry(pipe, view)
    d = np.load(reference_npz(pipe, geom['orientation']))
    h, z, Es, Ei, pos = lab_slice(d, plane, geom['orientation'])
    I0 = np.mean(sum(np.abs(Ei[c]) ** 2 for c in COMPS))
    I = {'total': sum(np.abs(Es[c] + Ei[c]) ** 2 for c in COMPS) / I0,
         'scattered': sum(np.abs(Es[c]) ** 2 for c in COMPS) / I0}
    n2d = phantom_slice(pipe, plane, h, z, pos)
    return geom, h, z, I, n2d, pos


def save_figure(fig, path_base, rv):
    if rv.save_png:
        fig.savefig(f'{path_base}.png', dpi=200)
    if rv.save_svg:
        fig.savefig(f'{path_base}.svg')


def save_lateral_reference_view(pipe, view):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rv = pipe.cfg.reference_views
    out = pipe.view_dir(view) / 'reference_views'
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for plane in rv.planes:
        geom, h, z, I, n2d, pos = case_images(pipe, view, plane)
        for mode in rv.field_modes:
            alpha = rv.arcsinh_alpha or auto_alpha(I[mode], n2d, z, pipe.n_medium, geom['camera'] == 'top')
            D = compute_display_intensity(I[mode], rv.intensity_mode, alpha)
            w = 6.4
            fig, ax = plt.subplots(figsize=(w, w * (z[-1] - z[0]) / (h[-1] - h[0]) + 1.6))
            title = (f"{TITLES.get(view, view)} — campo {'total' if mode == 'total' else 'dispersado'} — corte {plane}\n"
                     f"vista {view} · incidencia normal (kx = ky = 0) · λ = {pipe.wavelength:.4f} µm · "
                     f"n_medio = {pipe.n_medium:g} · {'x' if plane == 'YZ' else 'y'} = {pos:+.2f} µm\n"
                     f"{rv.intensity_mode} · {rv.colormap}")
            im = plot_lateral_field(ax, D, h, z, n2d, pipe.n_medium, geom, rv, title, alpha)
            cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
            cb.set_label(colorbar_label(rv, alpha), fontsize=7)
            cb.ax.tick_params(labelsize=7)
            fig.tight_layout()
            base = out / f'lateral_field_{mode}_{plane}_{rv.intensity_mode}_{rv.colormap}'
            save_figure(fig, base, rv)
            plt.close(fig)
            written.append(base)
    return written


def combined_figure(pipe, views):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rv = pipe.cfg.reference_views
    grid = [['top', 'bottom'], ['refl_top', 'refl_bottom']]
    if not all(v in views for row in grid for v in row):
        return []
    out = pipe.out / 'reference_views'
    written = []
    for plane in rv.planes:
        data = {v: case_images(pipe, v, plane) for row in grid for v in row}
        for mode in rv.field_modes:
            h0, z0 = data['top'][1], data['top'][2]
            w = 5.2
            fig, axs = plt.subplots(2, 2, figsize=(2 * w + 1.2, 2 * w * (z0[-1] - z0[0]) / (h0[-1] - h0[0]) + 1.8))
            for i, row in enumerate(grid):
                for j, v in enumerate(row):
                    geom, h, z, I, n2d, pos = data[v]
                    alpha = rv.arcsinh_alpha or auto_alpha(I[mode], n2d, z, pipe.n_medium, geom['camera'] == 'top')
                    D = compute_display_intensity(I[mode], rv.intensity_mode, alpha)
                    im = plot_lateral_field(axs[i, j], D, h, z, n2d, pipe.n_medium, geom, rv,
                                            f"{TITLES[v]} — corte {plane}" +
                                            ('' if rv.intensity_mode == 'raw' else f"  (α = {alpha:.2g})"), alpha)
            fig.suptitle(f"Campo {'total' if mode == 'total' else 'dispersado'} en el dominio FDTD, iluminación central "
                         f"(kx = ky = 0), λ = {pipe.wavelength:.4f} µm, n_medio = {pipe.n_medium:g}\n"
                         f"{rv.intensity_mode} · {rv.colormap} · contorno: phantom (índice simulado)", fontsize=10)
            cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02)
            cb.set_label(colorbar_label(rv, float('nan')).split(',  $\\alpha$')[0] + '  (α propio de cada panel)'
                         if rv.intensity_mode != 'raw' else colorbar_label(rv, 0), fontsize=8)
            base = out / f'lateral_reference_4casos_{mode}_{plane}_{rv.intensity_mode}_{rv.colormap}'
            save_figure(fig, base, rv)
            plt.close(fig)
            written.append(base)
    return written


def run_reference(pipe, views=None):
    views = list(views or (list(pipe.acq) + list(pipe.cfg.reflection.views)))
    for o in sorted({case_geometry(pipe, v)['orientation'] for v in views}):
        simulate_reference(pipe, o)
    written = []
    for v in views:
        written += save_lateral_reference_view(pipe, v)
    written += combined_figure(pipe, views)
    return written
