# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón

import math
import re
from pathlib import Path

import numpy as np
from PIL import Image


def read_config_manip(path):
    config = {}
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        parts = line.split('#', 1)[0].split()
        if len(parts) < 2:
            continue
        key, value = parts[:2]
        try:
            value = float(value)
        except ValueError:
            pass
        config[key.upper()] = value
    return config


def list_holograms(folder):
    files = [p for p in Path(folder).glob('*.pgm') if re.fullmatch(r'i\d+\.pgm', p.name)]
    return sorted(files, key=lambda p: int(p.stem[1:]))


def read_pgm(path):
    with Image.open(path) as img:
        return np.asarray(img, dtype=np.float64)


def tukey_window(n, r):
    if r <= 0:
        return np.ones(n)
    x = np.arange(n) / (n - 1)
    if r >= 1:
        return 0.5 * (1 - np.cos(2 * np.pi * x))
    w = np.ones(n)
    lo, hi = x < r / 2, x > 1 - r / 2
    w[lo] = 0.5 * (1 + np.cos(2 * np.pi / r * (x[lo] - r / 2)))
    w[hi] = 0.5 * (1 + np.cos(2 * np.pi / r * (x[hi] - 1 + r / 2)))
    return w


def coord_to_coord_shift(cx, cy, nx, ny):
    sx = cx - nx / 2 if cx - nx / 2 > 0 else nx / 2 + cx
    sy = cy - ny / 2 if cy - ny / 2 > 0 else ny / 2 + cy
    return int(math.floor(sx + 0.5)), int(math.floor(sy + 0.5))


def fourier_geometry(cfg, shape):
    ny, nx = shape
    gtot = (cfg['F_TUBE'] / cfg['F_OBJ']) / cfg['RF']
    dx = cfg['TPCAM'] / gtot
    fov = cfg['DIM_ROI'] * dx
    dk = 2 * np.pi / fov
    k0 = 2 * np.pi * cfg['N0'] / cfg['LAMBDA']
    rewald = k0 / dk
    fmax = math.ceil(rewald * cfg['NA'] / cfg['N0'])
    dim = 2 * fmax
    csx, csy = coord_to_coord_shift(cfg['CIRCLE_CX'], cfg['CIRCLE_CY'], nx, ny)
    rows = np.mod(np.arange(csy - dim // 2 + 1, csy + dim // 2 + 1) - 1, ny)
    cols = np.mod(np.arange(csx - dim // 2 + 1, csx + dim // 2 + 1) - 1, nx)
    return dict(Gtot=gtot, dx_obj=dx, FOV=fov, dk=dk, k0=k0, REwald=rewald,
                fmaxHolo=fmax, dimHolo=dim, CentreXShift=csx, CentreYShift=csy,
                rows=rows, cols=cols, patch_center_row=dim // 2 + 1,
                patch_center_col=dim // 2 + 1, camera_Nx=nx, camera_Ny=ny)


class Acquisition:
    def __init__(self, folder):
        self.folder = Path(folder)
        self.cfg = read_config_manip(self.folder / 'config_manip.txt')
        n = int(self.cfg['DIM_ROI'])
        self.shape = (n, n)
        self.geom = fourier_geometry(self.cfg, self.shape)
        self.wavelength_um = self.cfg['LAMBDA'] * 1e6
        self.n0 = self.cfg['N0']
        self.NA = self.cfg['NA']
        self.dx_cam_um = self.geom['dx_obj'] * 1e6
