# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón
# Objects. v1.2: medium index taken from the object file, no index correction, automatic crop.

import json
from pathlib import Path

import h5py
import numpy as np


def rotate_point(p, axis):
    x, y, z = p
    if axis == 'y':
        return (-x, y, -z)
    if axis == 'x':
        return (x, -y, -z)
    return (x, y, z)


def rotate_volume_zyx(n, axis):
    if axis == 'y':
        return n[::-1, :, ::-1]
    if axis == 'x':
        return n[::-1, ::-1, :]
    return n


def write_voxel_h5(path, n_zyx, voxel_um, metadata=None):
    n_zyx = np.asarray(n_zyx, dtype=np.float32)
    if n_zyx.ndim != 3:
        raise ValueError('Se espera un volumen n[z, y, x].')
    vox = np.broadcast_to(np.asarray(voxel_um, float), (3,))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, 'w') as f:
        f.create_dataset('n', data=n_zyx, compression='gzip')
        for axis, count, step in zip('xyz', n_zyx.shape[::-1], vox):
            f[f'coords/{axis}_um'] = (np.arange(count) - (count - 1) / 2) * step
        f.attrs['axis_order'] = 'zyx'
        f.attrs['metadata_json'] = json.dumps(metadata or {})


def read_voxel_h5(path):
    with h5py.File(path, 'r') as f:
        order = f.attrs.get('axis_order', 'zyx')
        if isinstance(order, bytes):
            order = order.decode()
        if not isinstance(order, str):
            order = str(np.asarray(order).ravel()[0])
        if order.strip().rstrip('\x00') != 'zyx':
            raise ValueError(f"axis_order debe ser 'zyx' (vista h5py), no {order!r}.")
        n = np.asarray(f['n'][()], dtype=float)
        coords = {a: np.asarray(f[f'coords/{a}_um'][()], dtype=float).ravel() for a in 'xyz'}
    if not np.all(np.isfinite(n)) or np.any(n <= 0):
        raise ValueError('n debe ser finito y positivo.')
    for axis, count in zip('zyx', n.shape):
        c = coords[axis]
        if len(c) != count:
            raise ValueError(f'coords/{axis}_um no coincide con el tamaño de n.')
        if count > 1 and not np.allclose(np.diff(c), c[1] - c[0], rtol=1e-6, atol=1e-9):
            raise ValueError(f'El muestreo en {axis} debe ser uniforme.')
    step = [float(coords[a][1] - coords[a][0]) if len(coords[a]) > 1 else 1.0 for a in 'xyz']
    return n, step


def _mat_value(group, key):
    ds = group[key]
    v = ds[()]
    if ds.attrs.get('MATLAB_class') == b'char':
        return ''.join(chr(int(c)) for c in v.flatten())
    return v.item() if v.size == 1 else v.flatten().tolist()


MAT_MEDIUM_KEYS = ('n_medium', 'phantom_n_immersion')


def medium_index_from_file(path):
    """Índice del medio que indica el propio objeto, o None.

    .mat (GT registrado): meta.n_medium, o si no existe meta.phantom_n_immersion.
    .h5 (voxel): atributo n_medium, o n_medium en metadata_json, o la mediana del borde.
    """
    path = Path(path)
    if not path.exists():
        return None
    with h5py.File(path, 'r') as f:
        if 'phantom_gt_dn' in f:
            for key in MAT_MEDIUM_KEYS:
                if key in f.get('meta', {}):
                    return float(_mat_value(f['meta'], key))
            return None
        if 'n_medium' in f.attrs:
            return float(np.asarray(f.attrs['n_medium']).ravel()[0])
        meta = json.loads(f.attrs.get('metadata_json', '{}') or '{}')
        if meta.get('n_medium'):
            return float(meta['n_medium'])
        n = f['n']
        border = np.concatenate([np.ravel(n[0]), np.ravel(n[-1]), np.ravel(n[:, 0]),
                                 np.ravel(n[:, -1]), np.ravel(n[:, :, 0]), np.ravel(n[:, :, -1])])
        return float(np.median(border))


def auto_crop(path, pad_um=0.8, threshold=None, chunk=32):
    """Caja del soporte de Δn (|Δn| > umbral) más pad_um por cada lado, en vóxeles del GT:
    [[x0, x1], [y0, y1], [z0, z1]], base 0, fin excluido. Se lee por bloques en z."""
    with h5py.File(path, 'r') as f:
        ds = f['phantom_gt_dn']
        meta = f.get('meta', {})
        if threshold is None:
            threshold = float(_mat_value(meta, 'mask_threshold')) if 'mask_threshold' in meta else 1e-4
        voxel = np.asarray(_mat_value(meta, 'rec_voxel_um_xyz'), dtype=float)
        nz, nx, ny = ds.shape
        lo, hi = np.array([nz, nx, ny]), np.zeros(3, int)
        for z in range(0, nz, chunk):
            idx = np.argwhere(np.abs(ds[z:z + chunk]) > threshold)
            if len(idx):
                idx[:, 0] += z
                lo, hi = np.minimum(lo, idx.min(0)), np.maximum(hi, idx.max(0) + 1)
    if np.any(hi <= lo):
        raise ValueError(f'{path}: no hay vóxeles con |Δn| > {threshold}.')
    pad = np.ceil(pad_um / np.broadcast_to(voxel, (3,))[[2, 0, 1]]).astype(int)
    lo, hi = np.maximum(lo - pad, 0), np.minimum(hi + pad, [nz, nx, ny])
    (z0, x0, y0), (z1, x1, y1) = lo.tolist(), hi.tolist()
    return [[x0, x1], [y0, y1], [z0, z1]]


def registered_gt_crop(path, crop_xyz=None, n_medium=None, n_immersion_gt=None, soft_mask_dn=0.01,
                       pad_um=0.8):
    """GT registrado (Δn) recortado: n = n_medio + Δn, sin corrección.

    n_medio: el indicado (config) o, si es None, el del propio archivo (medium_index_from_file).
    crop_xyz None o 'auto': caja del soporte + pad_um. n_immersion_gt solo se aplica si se pide
    explícitamente (compatibilidad v1.1: Δn referido a otra inmersión, máscara suave).
    """
    if n_medium is None:
        n_medium = medium_index_from_file(path)
        if n_medium is None:
            raise ValueError(f'{path}: el archivo no indica el índice del medio (meta.n_medium); '
                             'indícalo en object.n_medium.')
    if crop_xyz in (None, 'auto', []):
        crop_xyz = auto_crop(path, pad_um)
    with h5py.File(path, 'r') as f:
        ds = f['phantom_gt_dn']
        voxel = np.asarray(_mat_value(f['meta'], 'rec_voxel_um_xyz'), dtype=float)
        (x0, x1), (y0, y1), (z0, z1) = crop_xyz
        dn = np.ascontiguousarray(ds[z0:z1, x0:x1, y0:y1].transpose(0, 2, 1), dtype=np.float64)
        nz, nx, ny = ds.shape
    n = n_medium + dn
    if n_immersion_gt is not None:
        n += (n_immersion_gt - n_medium) * np.clip(dn / soft_mask_dn, 0, 1)
    border = max(np.abs(dn[[0, -1]]).max(), np.abs(dn[:, [0, -1]]).max(), np.abs(dn[:, :, [0, -1]]).max())
    shape_xyz = np.array([nx, ny, nz])
    bounds = np.array(crop_xyz)
    center = ((bounds[:, 0] + bounds[:, 1] - 1) / 2 - (shape_xyz - 1) / 2) * voxel
    meta = dict(source=str(Path(path).resolve()), crop_xyz=bounds.tolist(), n_medium=n_medium,
                n_immersion_gt=n_immersion_gt, soft_mask_dn=soft_mask_dn,
                crop_center_um_xyz=center.tolist(), border_max_dn=float(border),
                n_min=float(n.min()), n_max=float(n.max()))
    return n.astype(np.float32), voxel, meta


def sphere_volume(radius_um, dn, n_medium, voxel_um=0.1, box_um=None, edge_um=0.05):
    box = box_um or 2 * radius_um + 1.0
    m = int(round(box / voxel_um))
    g = (np.arange(m) - (m - 1) / 2) * voxel_um
    Z, Y, X = np.meshgrid(g, g, g, indexing='ij')
    r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    return (n_medium + dn * 0.5 * (1 - np.tanh((r - radius_um) / edge_um))).astype(np.float32)
