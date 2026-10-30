# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# GUI data access: objects, illumination lists, holograms, sensor fields, saved figures (no Qt).

import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np

from fdtd_meep.acquisition import Acquisition, list_holograms, read_config_manip, read_pgm
from fdtd_meep.config import ROOT
from fdtd_meep.illuminations import read_illumination_list
from fdtd_meep.objects import (auto_crop, medium_index_from_file, read_voxel_h5, registered_gt_crop,
                               sphere_volume)

TRANSMISSION_VIEWS = ('bottom', 'top')
REFLECTION_VIEWS = ('refl_top', 'refl_bottom')
FREE_VIEWS = ('transmission', 'reflection')


def resolve(path):
    if not path:
        return None
    p = Path(path).expanduser()
    return p if p.is_absolute() else ROOT / p


def relative(path):
    """Ruta relativa a la raíz del proyecto si está dentro; si no, absoluta."""
    p = Path(path).expanduser().resolve()
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def output_dir(cfg):
    return resolve(cfg.get('output') or f"outputs/{cfg.get('name', 'run')}")


def is_free(cfg):
    return not cfg.get('views')


def acquisition_info(folder):
    """Resumen de una carpeta de adquisición (config_manip.txt + iNNN.pgm) o None."""
    folder = resolve(folder)
    if not folder or not (folder / 'config_manip.txt').exists():
        return None
    c = read_config_manip(folder / 'config_manip.txt')
    return dict(wavelength_um=c.get('LAMBDA', float('nan')) * 1e6, n0=c.get('N0'), NA=c.get('NA'),
                dim=int(c.get('DIM_ROI', 0)), n_holograms=len(list_holograms(folder)))


def medium_info(cfg):
    """(n_medio, origen) con la misma prioridad que el pipeline: config → objeto → N0."""
    o = cfg.get('object', {})
    if o.get('n_medium'):
        return float(o['n_medium']), 'set in config'
    if o.get('source') in ('registered_gt', 'voxel_h5') and resolve(o.get('file')):
        try:
            n = medium_index_from_file(resolve(o['file']))
        except (OSError, KeyError, ValueError):
            n = None
        if n:
            return n, 'from object file'
    for folder in cfg.get('views', {}).values():
        info = acquisition_info(folder)
        if info:
            return float(info['n0']), 'from config_manip N0'
    return None, 'not set'


def n_medium_of(cfg):
    return medium_info(cfg)[0]


@lru_cache(maxsize=8)
def _auto_crop(path, mtime, pad_um):
    return auto_crop(path, pad_um)


def crop_for(cfg):
    """crop_xyz de la configuración o, si está vacío, el automático (soporte + margen)."""
    o = cfg.get('object', {})
    if o.get('crop_xyz') and o['crop_xyz'] != 'auto':
        return o['crop_xyz']
    p = resolve(o.get('file'))
    return _auto_crop(str(p), p.stat().st_mtime, float(o.get('crop_pad_um', 0.8)))


def load_object(cfg):
    """Volumen n[z, y, x] y paso (x, y, z) en µm: el ya preparado, si existe, o el de entrada."""
    prepared = output_dir(cfg) / 'object' / 'object.h5'
    if prepared.exists():
        n, step = read_voxel_h5(prepared)
        return n, step, f'prepared ({prepared.parents[1].name}/object/object.h5)'
    o = cfg.get('object', {})
    source = o.get('source', 'registered_gt')
    if source == 'sphere':
        nm = n_medium_of(cfg) or 1.5
        vox = float(o.get('voxel_um', 0.1))
        n = sphere_volume(float(o.get('radius_um', 1.5)), float(o.get('dn', 0.03)), nm, vox,
                          edge_um=float(o.get('edge_um', 0.05)))
        return n, [vox] * 3, 'sphere preview'
    path = resolve(o.get('file'))
    if not path or not path.exists():
        raise FileNotFoundError(f'Object file not found: {o.get("file")}')
    if source == 'voxel_h5':
        n, step = read_voxel_h5(path)
        return n, step, path.name
    n, vox, _ = registered_gt_crop(path, crop_for(cfg), n_medium_of(cfg), o.get('n_immersion_gt'),
                                   o.get('soft_mask_dn', 0.01))
    return n, list(np.broadcast_to(vox, (3,))), f'{path.name} (crop)'


def object_slice(n, step, plane):
    """Corte central: XZ, YZ o XY. Devuelve (imagen, extent, xlabel, ylabel)."""
    nz, ny, nx = n.shape
    ext = [((m - 1) / 2 + 0.5) * s for m, s in zip((nx, ny, nz), step)]
    if plane == 'XZ':
        return n[:, ny // 2, :], [-ext[0], ext[0], -ext[2], ext[2]], 'x (µm)', 'z (µm)'
    if plane == 'YZ':
        return n[:, :, nx // 2], [-ext[1], ext[1], -ext[2], ext[2]], 'y (µm)', 'z (µm)'
    return n[nz // 2], [-ext[0], ext[0], -ext[1], ext[1]], 'x (µm)', 'y (µm)'


def _read_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def illumination_sets(cfg):
    """Conjuntos a dibujar en (kx, ky): [(etiqueta, kx, ky, nombres, tipo)], tipo 'all' o 'sel'."""
    out = []
    if is_free(cfg):
        p = resolve(cfg.get('illuminations'))
        if p and p.exists():
            rows = read_illumination_list(p)
            out.append(('list', [r['kx'] for r in rows], [r['ky'] for r in rows],
                        [r['name'] for r in rows], 'sel'))
        return out
    lists = cfg.get('selection', {}).get('lists', {})
    for view in cfg.get('views', {}):
        sel = output_dir(cfg) / view / 'selection'
        if (sel / 'illumination_positions_all.csv').exists():
            rows = _read_csv(sel / 'illumination_positions_all.csv')
            out.append((f'{view}: all', [float(r['kx_norm_raw']) for r in rows],
                        [float(r['ky_norm_raw']) for r in rows], [r['filename'] for r in rows], 'all'))
        if (sel / 'illumination_positions.csv').exists():
            rows = _read_csv(sel / 'illumination_positions.csv')
            out.append((f'{view}: selected', [float(r['kx_norm_raw']) for r in rows],
                        [float(r['ky_norm_raw']) for r in rows], [r['new_filename'][:-4] for r in rows], 'sel'))
        elif lists.get(view) and resolve(lists[view]).exists():
            rows = [r for r in read_illumination_list(resolve(lists[view])) if 'kx' in r]
            if rows:
                out.append((f'{view}: list', [r['kx'] for r in rows], [r['ky'] for r in rows],
                            [r['name'] for r in rows], 'sel'))
    return out


def available_views(cfg):
    out_dir = output_dir(cfg)
    if is_free(cfg):
        names = list(FREE_VIEWS)
    else:
        names = list(cfg.get('views', {})) + list(cfg.get('reflection', {}).get('views', REFLECTION_VIEWS))
    return [v for v in names if (out_dir / v).exists() or v in cfg.get('views', {})]


def hologram_sources(cfg, view):
    """{etiqueta: carpeta} con hologramas iNNN.pgm para esta vista."""
    vdir = output_dir(cfg) / view
    out = {}
    if view in cfg.get('views', {}):
        folder = resolve(cfg['views'][view])
        if folder and folder.exists():
            out['input'] = folder
    for d in sorted(vdir.glob('real_subset')) + sorted(vdir.glob('sim_subset*')):
        if d.is_dir() and list_holograms(d):
            out[d.name.replace('_subset', '')] = d
    return out


def illumination_names(cfg, view, folder=None):
    names = set()
    if folder:
        names |= {p.stem for p in list_holograms(folder)}
    fdir = output_dir(cfg) / view / 'fields'
    if fdir.exists():
        names |= {p.stem for p in fdir.glob('*.npz')}
    return sorted(names)


def hologram_maps(folder, name, tukey=0.1, annulus=(17.0, 33.0), order=4):
    """Intensidad del holograma y campo demodulado del orden +1 (fondo de fase corregido)."""
    from skimage.restoration import unwrap_phase
    from fdtd_meep.camera import demodulate, remove_background
    from fdtd_meep.specular import analyze_hologram, window2d
    folder = Path(folder)
    acq = Acquisition(folder)
    holo = read_pgm(folder / f'{name}.pgm')
    pk = analyze_hologram(holo, acq.geom, window2d(acq.shape[1], tukey))
    U, dx = demodulate(holo, acq.geom, pk['drow_subpix'], pk['dcol_subpix'], tukey)
    try:
        U = remove_background(U, dx, annulus, order)
    except Exception:
        pass
    phase = np.angle(U)
    h_cam = acq.shape[0] * acq.dx_cam_um / 2
    h_pat = U.shape[0] * dx / 2
    return dict(intensity=(holo, [-h_cam, h_cam, -h_cam, h_cam]),
                phase=(phase, [-h_pat, h_pat, -h_pat, h_pat]),
                unwrapped=(unwrap_phase(phase), [-h_pat, h_pat, -h_pat, h_pat]),
                amplitude=(np.abs(U), [-h_pat, h_pat, -h_pat, h_pat]), peak=pk)


def field_maps(npz_path):
    """Campo relativo en el sensor: |1 + S| y fase en transmisión, |S| y fase en reflexión."""
    d = np.load(npz_path)
    pol = np.abs(d['pol']) if 'pol' in d else np.array([1, 0, 0])
    comp = 'Ex' if pol[0] >= pol[1] else 'Ey'
    S = d[f'Es_{comp}'] / d[f'Einc_{comp}']
    reflection = bool(d['reflection']) if 'reflection' in d else False
    U = S if reflection else 1 + S
    x, y = d['x'], d['y']
    return dict(amplitude=np.abs(U), phase=np.angle(U), extent=[x[0], x[-1], y[0], y[-1]],
                label='|S|' if reflection else '|1 + S|', reflection=reflection, component=comp,
                z=float(d['z']), kx_norm=float(d['kx_norm']) if 'kx_norm' in d else float('nan'),
                ky_norm=float(d['ky_norm']) if 'ky_norm' in d else float('nan'),
                n_cells_M=float(d['n_cells']) / 1e6 if 'n_cells' in d else float('nan'))


def report_rows(cfg, view):
    """Filas del informe de la vista (sim_report.csv o free_report.csv), por nombre sin .pgm."""
    vdir = output_dir(cfg) / view
    for p in (vdir / 'sim_subset' / 'sim_report.csv', vdir / 'free_report.csv'):
        if p.exists():
            rows = _read_csv(p)
            return {(r.get('new_filename') or r.get('name', '')).replace('.pgm', ''): r for r in rows}
    sel = vdir / 'selection' / 'illumination_positions.csv'
    if sel.exists():
        return {r['new_filename'][:-4]: r for r in _read_csv(sel)}
    return {}


def calibration_summary(cfg):
    p = output_dir(cfg) / 'calibration' / 'calibration.json'
    if not p.exists():
        return None
    cal = json.loads(p.read_text())
    return {v: cal[v]['best'] for v in cal if isinstance(cal[v], dict) and 'best' in cal[v]}


def saved_figures(cfg):
    out = output_dir(cfg)
    if not out.exists():
        return []
    return sorted(p for p in out.rglob('*.png') if 'fields' not in p.parts)
