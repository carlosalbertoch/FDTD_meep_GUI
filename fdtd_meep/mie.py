# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón

import numpy as np


def _rotation_to_z(khat):
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(khat, z)
    s, c = np.linalg.norm(v), float(np.dot(khat, z))
    if s < 1e-12:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / s ** 2)


def angular_spectrum(kx, ky, k_inc, pol, n_sphere, n_medium, radius_um, wavelength_um, backward=False):
    import miepython
    k = 2 * np.pi * n_medium / wavelength_um
    kz2 = k ** 2 - kx ** 2 - ky ** 2
    prop = kz2 > 0
    kz = np.sqrt(np.where(prop, kz2, 1.0))
    ks = np.stack([kx, ky, -kz if backward else kz], -1) / k
    ki = np.asarray(k_inc, float) / np.linalg.norm(k_inc)
    R = _rotation_to_z(ki)
    e = R @ np.asarray(pol, complex)
    s = ks @ R.T
    theta = np.arccos(np.clip(s[..., 2], -1, 1))
    phi = np.arctan2(s[..., 1], s[..., 0])
    S1, S2 = miepython.S1_S2(n_sphere / n_medium, k * radius_um, np.cos(theta).ravel(), norm='wiscombe')
    S1, S2 = np.conj(S1.reshape(theta.shape)), np.conj(S2.reshape(theta.shape))
    cp, sp, ct, st = np.cos(phi), np.sin(phi), np.cos(theta), np.sin(theta)
    E_par = e[0] * cp + e[1] * sp
    E_perp = e[0] * sp - e[1] * cp
    e_theta = np.stack([ct * cp, ct * sp, -st], -1)
    e_phi = np.stack([-sp, cp, np.zeros_like(cp)], -1)
    f_rot = (S2 * E_par)[..., None] * e_theta - (S1 * E_perp)[..., None] * e_phi
    f = (f_rot @ R) / (-1j * k)
    return np.where(prop[..., None], 2j * np.pi * f / kz[..., None], 0)
