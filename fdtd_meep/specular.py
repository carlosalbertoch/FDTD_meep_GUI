# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón
# Specular peaks and selection. v1.2: selection can follow an illumination list.

import csv
import math
import re
import shutil
from pathlib import Path

import numpy as np

from .acquisition import Acquisition, list_holograms, read_pgm, tukey_window

CENTROID_RADIUS = 3
COLUMNS = ['index', 'filename', 'peak_row', 'peak_col', 'peak_row_subpix', 'peak_col_subpix',
           'drow_subpix', 'dcol_subpix', 'peak_quality', 'peak_snr', 'kx_norm_raw', 'ky_norm_raw',
           'rho_norm', 'theta_deg', 'phi_deg', 'fps_rank']


def window2d(n, tukey):
    w = tukey_window(n, tukey)
    return np.sqrt(np.outer(w, w))


def patch_power(holo, geom, window):
    F = np.fft.fft2(window * holo)
    return np.abs(F[np.ix_(geom['rows'], geom['cols'])]) ** 2


def analyze_hologram(holo, geom, window):
    P = patch_power(holo, geom, window)
    pr, pc = np.unravel_index(int(np.argmax(P)), P.shape)
    nr, nc = P.shape
    r1, r2 = max(0, pr - CENTROID_RADIUS), min(nr - 1, pr + CENTROID_RADIUS)
    c1, c2 = max(0, pc - CENTROID_RADIUS), min(nc - 1, pc + CENTROID_RADIUS)
    loc = P[r1:r2 + 1, c1:c2 + 1]
    R, C = np.mgrid[r1:r2 + 1, c1:c2 + 1]
    border = np.zeros(loc.shape, bool)
    border[[0, -1], :] = True
    border[:, [0, -1]] = True
    bg = float(np.median(loc[border]))
    W = np.maximum(loc - bg, 0)
    sw = W.sum()
    rc, cc = (pr, pc) if sw <= np.finfo(float).eps else ((R * W).sum() / sw, (C * W).sum() / sw)
    eps = np.finfo(float).eps
    center = geom['patch_center_row']
    drow, dcol = rc + 1 - center, cc + 1 - center
    kx, ky = dcol / geom['REwald'], drow / geom['REwald']
    rho = math.hypot(kx, ky)
    return dict(peak_row=pr + 1, peak_col=pc + 1, peak_row_subpix=rc + 1, peak_col_subpix=cc + 1,
                drow_subpix=drow, dcol_subpix=dcol,
                peak_quality=float(P[pr, pc] / (P.mean() + eps)),
                peak_snr=float(P[pr, pc] / (bg + eps)),
                kx_norm_raw=kx, ky_norm_raw=ky, rho_norm=rho,
                theta_deg=math.degrees(math.asin(min(rho, 1))),
                phi_deg=math.degrees(math.atan2(ky, kx)))


def fps_order(kx, ky):
    K = np.column_stack([kx, ky])
    n = len(K)
    order = [int(np.argmin(np.hypot(kx, ky)))]
    dmin = np.hypot(*(K - K[order[0]]).T)
    for _ in range(1, n):
        dmin[order] = -np.inf
        j = int(np.argmax(dmin))
        order.append(j)
        dmin = np.minimum(dmin, np.hypot(*(K - K[j]).T))
    return order


def extract_all(folder, tukey):
    acq = Acquisition(folder)
    window = window2d(acq.shape[1], tukey)
    rows = []
    for path in list_holograms(folder):
        r = analyze_hologram(read_pgm(path), acq.geom, window)
        rows.append(dict(index=int(path.stem[1:]), filename=path.name, **r))
    order = fps_order(np.array([r['kx_norm_raw'] for r in rows]),
                      np.array([r['ky_norm_raw'] for r in rows]))
    for rank, j in enumerate(order, 1):
        rows[j]['fps_rank'] = rank
    return rows, order


def write_csv(path, rows, columns):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def read_selection(view_dir):
    with open(Path(view_dir) / 'selection' / 'illumination_positions.csv') as f:
        return list(csv.DictReader(f))


def select_view(view, folder, out_view_dir, n_sel, tukey, wanted=None):
    """wanted: lista de iluminaciones (read_illumination_list); sustituye a la selección FPS."""
    from scipy.io import savemat
    from .illuminations import match_to_holograms
    sel_dir = Path(out_view_dir) / 'selection'
    sel_dir.mkdir(parents=True, exist_ok=True)
    rows, order = extract_all(folder, tukey)
    write_csv(sel_dir / 'illumination_positions_all.csv', rows, COLUMNS)
    picks = order[:n_sel] if wanted is None else match_to_holograms(rows, wanted)
    n_sel = len(picks)
    how = 'orden FPS en (kx,ky) RAW' if wanted is None else 'lista de iluminaciones'
    chosen = []
    for new_index, j in enumerate(picks):
        chosen.append(dict(rows[j], new_index=new_index, new_filename=f'i{new_index:03d}.pgm'))
    cols = COLUMNS + ['new_index', 'new_filename']
    write_csv(sel_dir / 'illumination_positions.csv', chosen, cols)
    savemat(sel_dir / 'illumination_positions.mat',
            {'illum': {k: np.array([r[k] for r in chosen], dtype=object if 'filename' in k else None)
                       for k in cols},
             'meta': {'view': view, 'folder': str(folder), 'n_selected': n_sel,
                      'n_total': len(rows), 'tukey': tukey,
                      'k_axis_convention': 'RAW: +dcol -> +kx, +drow -> +ky'}})
    with open(sel_dir / 'seleccion.txt', 'w') as f:
        f.write(f'# {Path(folder).name}: {n_sel} de {len(rows)} hologramas, {how}\n')
        f.write('# nuevo_nombre  archivo_original  kx_norm_raw  ky_norm_raw  theta_deg\n')
        for r in chosen:
            f.write(f"{r['new_filename']}  {r['filename']}  {r['kx_norm_raw']:+.5f}  "
                    f"{r['ky_norm_raw']:+.5f}  {r['theta_deg']:6.2f}\n")
    real = Path(out_view_dir) / 'real_subset'
    if real.exists():
        shutil.rmtree(real)
    real.mkdir(parents=True)
    for r in chosen:
        shutil.copy(Path(folder) / r['filename'], real / r['new_filename'])
    for name in ('Intensite_ref.pgm', 'recon.txt'):
        if (Path(folder) / name).exists():
            shutil.copy(Path(folder) / name, real / name)
    txt = (Path(folder) / 'config_manip.txt').read_text()
    (real / 'config_manip.txt').write_text(re.sub(r'(?m)^NB_HOLO\s+\d+', f'NB_HOLO {n_sel}', txt))
    shutil.copy(sel_dir / 'seleccion.txt', real / 'seleccion.txt')
    coverage_figure(sel_dir / 'seleccion_k.png', rows, chosen, Path(folder).name,
                    'FPS' if wanted is None else 'lista')
    return chosen


def coverage_figure(path, rows, chosen, title, how='FPS'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot([r['kx_norm_raw'] for r in rows], [r['ky_norm_raw'] for r in rows], '.', color='0.7')
    ax.plot([r['kx_norm_raw'] for r in chosen], [r['ky_norm_raw'] for r in chosen], 'ro', ms=4)
    for r in chosen:
        ax.annotate(r['new_filename'][:-4], (r['kx_norm_raw'], r['ky_norm_raw']), fontsize=6)
    ax.set_aspect('equal')
    ax.grid(alpha=0.3)
    ax.set_xlabel('kx/k (RAW)')
    ax.set_ylabel('ky/k (RAW)')
    ax.set_title(f'{title}: {len(chosen)} de {len(rows)} iluminaciones ({how})')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
