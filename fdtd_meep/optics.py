# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón

import math

import numpy as np

from .acquisition import tukey_window


def tukey2d(ny, nx, alpha=0.25):
    return np.outer(tukey_window(ny, alpha), tukey_window(nx, alpha))


def filter_and_propagate(U, dx, wavelength_um, n_medium, NA=None, dz=0.0, pad=2):
    ny, nx = U.shape
    Py, Px = int(ny * pad), int(nx * pad)
    buf = np.zeros((Py, Px), dtype=complex)
    buf[:ny, :nx] = U
    KX, KY = np.meshgrid(2 * np.pi * np.fft.fftfreq(Px, dx), 2 * np.pi * np.fft.fftfreq(Py, dx))
    spec = np.fft.fft2(buf)
    if NA is not None:
        spec *= (KX ** 2 + KY ** 2) <= (2 * np.pi * NA / wavelength_um) ** 2
    if dz:
        k = 2 * np.pi * n_medium / wavelength_um
        kz2 = k ** 2 - KX ** 2 - KY ** 2
        prop = kz2 > 0
        spec *= np.where(prop, np.exp(1j * np.sqrt(np.where(prop, kz2, 0)) * dz), 0)
    return np.fft.ifft2(spec)[:ny, :nx]


def focus_field(S_rel, x, y, kx, ky, dz, NA, wavelength_um, n_medium, taper=0.15):
    X, Y = np.meshgrid(x, y)
    ramp = np.exp(1j * (kx * X + ky * Y))
    kz = math.sqrt((2 * np.pi * n_medium / wavelength_um) ** 2 - kx ** 2 - ky ** 2)
    S = S_rel * tukey2d(*S_rel.shape, alpha=taper) * ramp
    S = filter_and_propagate(S, x[1] - x[0], wavelength_um, n_medium, NA=NA, dz=dz)
    return S / (ramp * np.exp(1j * kz * dz))
