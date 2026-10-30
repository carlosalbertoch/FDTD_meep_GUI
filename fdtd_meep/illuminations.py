# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón
# Illumination list files (CSV): read, generate and match to real holograms.

import csv
import math
from pathlib import Path

import numpy as np

ALIASES = {'kx_norm_raw': 'kx', 'ky_norm_raw': 'ky', 'kx_norm': 'kx', 'ky_norm': 'ky',
           'theta': 'theta_deg', 'phi': 'phi_deg', 'new_filename': 'name'}


def read_illumination_list(path):
    """Lee una lista de iluminaciones. Columnas (cabecera obligatoria, '#' = comentario):

    - kx, ky: cosenos directores transversales kx/k, ky/k (k = 2π·n_medio/λ), o bien
    - theta_deg, phi_deg: kx = sin θ cos φ, ky = sin θ sin φ, o bien
    - filename: holograma real (iNNN.pgm), solo con adquisiciones reales.
    - name (opcional): nombre de salida; por defecto i000, i001… en el orden del archivo.

    illumination_positions.csv (salida de select) también se acepta tal cual.
    """
    with open(path, newline='') as f:
        lines = [l for l in f if l.strip() and not l.lstrip().startswith('#')]
    reader = csv.DictReader(lines)
    reader.fieldnames = [ALIASES.get(c.strip().lower(), c.strip().lower()) for c in reader.fieldnames]
    out = []
    for i, r in enumerate(reader):
        r = {k: (v or '').strip() for k, v in r.items() if k}
        row = {'name': r.get('name') or f'i{i:03d}'}
        row['name'] = row['name'][:-4] if row['name'].endswith('.pgm') else row['name']
        if r.get('filename'):
            row['filename'] = r['filename']
        if r.get('kx') and r.get('ky'):
            row['kx'], row['ky'] = float(r['kx']), float(r['ky'])
        elif r.get('theta_deg') and r.get('phi_deg'):
            t, p = math.radians(float(r['theta_deg'])), math.radians(float(r['phi_deg']))
            row['kx'], row['ky'] = math.sin(t) * math.cos(p), math.sin(t) * math.sin(p)
        elif 'filename' not in row:
            raise ValueError(f'{path}, fila {i + 2}: faltan kx,ky, theta_deg,phi_deg o filename.')
        if 'kx' in row:
            rho = math.hypot(row['kx'], row['ky'])
            if rho >= 1:
                raise ValueError(f'{path}, fila {i + 2}: kx² + ky² debe ser < 1 (vale {rho ** 2:.3f}).')
            row['theta_deg'] = math.degrees(math.asin(rho))
            row['phi_deg'] = math.degrees(math.atan2(row['ky'], row['kx']))
        out.append(row)
    if not out:
        raise ValueError(f'{path}: la lista está vacía.')
    names = [r['name'] for r in out]
    if len(set(names)) != len(names):
        raise ValueError(f'{path}: nombres repetidos en la columna name.')
    return out


def write_illumination_list(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='') as f:
        f.write('# kx, ky = kx/k, ky/k (k = 2*pi*n_medium/lambda); theta/phi in degrees, informative\n')
        w = csv.writer(f)
        w.writerow(['name', 'kx', 'ky', 'theta_deg', 'phi_deg'])
        for r in rows:
            w.writerow([r['name'], f"{r['kx']:.6f}", f"{r['ky']:.6f}",
                        f"{r['theta_deg']:.3f}", f"{r['phi_deg']:.3f}"])


def generate_illuminations(n, theta_max_deg, pattern='spiral'):
    """n direcciones dentro del cono θ ≤ θ_max, empezando por la normal.

    spiral: espiral de Fibonacci uniforme en (kx, ky); rings: anillos concéntricos.
    """
    rho_max = math.sin(math.radians(theta_max_deg))
    if pattern == 'spiral':
        i = np.arange(n)
        rho = rho_max * np.sqrt(i / max(n - 1, 1))
        phi = i * math.pi * (3 - math.sqrt(5))
    elif pattern == 'rings':
        n_rings = max(1, int(round(math.sqrt(n / math.pi))))
        weights = np.arange(1, n_rings + 1, dtype=float)
        per = np.maximum(1, np.round((n - 1) * weights / weights.sum())).astype(int)
        per[-1] += n - 1 - per.sum()
        rho, phi = [0.0], [0.0]
        for j, m in enumerate(per, 1):
            rho += [rho_max * j / n_rings] * m
            phi += list(2 * np.pi * np.arange(m) / m + (j % 2) * np.pi / m)
        rho, phi = np.array(rho), np.array(phi)
    else:
        raise ValueError(f'Patrón desconocido: {pattern}')
    rows = []
    for j, (r, p) in enumerate(zip(rho, phi)):
        kx, ky = float(r * math.cos(p)), float(r * math.sin(p))
        rows.append(dict(name=f'i{j:03d}', kx=kx, ky=ky, theta_deg=math.degrees(math.asin(r)),
                         phi_deg=math.degrees(math.atan2(ky, kx))))
    return rows


def match_to_holograms(rows, wanted):
    """Índices en rows (picos de los hologramas reales) de cada iluminación de wanted: por
    nombre de archivo, o por el pico especular más cercano en (kx, ky) RAW, sin repetir."""
    by_name = {r['filename']: i for i, r in enumerate(rows)}
    K = np.array([[r['kx_norm_raw'], r['ky_norm_raw']] for r in rows])
    picks = []
    for w in wanted:
        if 'filename' in w:
            if w['filename'] not in by_name:
                raise ValueError(f"{w['filename']} no está entre los hologramas de la vista.")
            j = by_name[w['filename']]
        else:
            d = np.hypot(K[:, 0] - w['kx'], K[:, 1] - w['ky'])
            d[picks] = np.inf
            j = int(np.argmin(d))
        if j in picks:
            raise ValueError(f"{rows[j]['filename']} aparece dos veces en la lista.")
        picks.append(j)
    return picks
