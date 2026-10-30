# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón

import numpy as np

from .specular import window2d

DIHEDRAL = {
    'id': ((1, 0), (0, 1)), 'flip_x': ((-1, 0), (0, 1)), 'flip_y': ((1, 0), (0, -1)),
    'rot180': ((-1, 0), (0, -1)), 'transpose': ((0, 1), (1, 0)),
    'rot90': ((0, -1), (1, 0)), 'rot270': ((0, 1), (-1, 0)), 'antitranspose': ((0, -1), (-1, 0)),
}


def to_camera_grid(S, x, y, M, t_um, dim, dx, offset=1.0):
    from scipy.interpolate import RegularGridInterpolator
    c = (np.arange(dim) - dim // 2) * dx
    Xc, Yc = np.meshgrid(c, c)
    Mi = np.array(M, float).T
    xo = Mi[0, 0] * (Xc - t_um[0]) + Mi[0, 1] * (Yc - t_um[1])
    yo = Mi[1, 0] * (Xc - t_um[0]) + Mi[1, 1] * (Yc - t_um[1])
    pts = np.stack([yo.ravel(), xo.ravel()], -1)
    fr = RegularGridInterpolator((y, x), S.real, bounds_error=False, fill_value=0.0)
    fi = RegularGridInterpolator((y, x), S.imag, bounds_error=False, fill_value=0.0)
    return offset + (fr(pts) + 1j * fi(pts)).reshape(Xc.shape)


def demodulate(holo, geom, drow, dcol, tukey=0.10):
    P = np.fft.fft2(window2d(holo.shape[1], tukey) * holo)[np.ix_(geom['rows'], geom['cols'])]
    dim = P.shape[0]
    c0 = geom['patch_center_row'] - 1
    ir, ic = int(round(drow)), int(round(dcol))
    P = np.roll(P, (-ir, -ic), (0, 1))
    Y, X = np.mgrid[:dim, :dim] - c0
    P = P * (np.hypot(X, Y) <= geom['fmaxHolo'])
    U = np.fft.ifft2(np.fft.ifftshift(P))
    Yn, Xn = np.mgrid[:dim, :dim]
    U = U * np.exp(-2j * np.pi * ((drow - ir) * Yn + (dcol - ic) * Xn) / dim)
    return U, geom['FOV'] * 1e6 / dim


def _poly_basis(n, order):
    y, x = (np.mgrid[:n, :n] - n / 2) / (n / 2)
    return np.stack([x ** i * y ** j for i in range(order + 1) for j in range(order + 1 - i)], -1)


def remove_background(U, dx, annulus_um=(17.0, 33.0), order=4, iterations=3):
    from skimage.restoration import unwrap_phase
    n = U.shape[0]
    phi = unwrap_phase(np.angle(U))
    la = np.log(np.abs(U) + 1e-12)
    yy, xx = (np.mgrid[:n, :n] - n // 2) * dx
    r = np.hypot(xx, yy)
    use = (r > annulus_um[0]) & (r < annulus_um[1])
    B = _poly_basis(n, order).reshape(n * n, -1)
    for _ in range(iterations):
        cp, *_ = np.linalg.lstsq(B[use.ravel()], phi[use], rcond=None)
        ca, *_ = np.linalg.lstsq(B[use.ravel()], la[use], rcond=None)
        rp = phi - (B @ cp).reshape(n, n)
        ra = la - (B @ ca).reshape(n, n)
        use &= (np.abs(rp) < 3 * rp[use].std()) & (np.abs(ra) < 3 * ra[use].std())
    rp -= np.median(rp[use])
    return np.exp(ra - np.median(ra[use]) + 1j * rp)


def match(T_real, S, x, y, dim, dx, transforms=None):
    A = T_real - 1
    fA = np.fft.fft2(A)
    out = []
    for name in (transforms or DIHEDRAL):
        M = DIHEDRAL[name]
        B0 = to_camera_grid(S, x, y, M, (0.0, 0.0), dim, dx) - 1
        for conj in (False, True):
            B = np.conj(B0) if conj else B0
            C = np.fft.ifft2(fA * np.conj(np.fft.fft2(B)))
            i = np.unravel_index(int(np.argmax(np.abs(C))), C.shape)
            sh = [(v if v < s / 2 else v - s) for v, s in zip(i, C.shape)]
            t = (sh[1] * dx, sh[0] * dx)
            Bt = to_camera_grid(S, x, y, M, t, dim, dx) - 1
            Bt = np.conj(Bt) if conj else Bt
            corr = abs(np.vdot(Bt, A)) / (np.linalg.norm(A) * np.linalg.norm(Bt))
            out.append(dict(M=name, conj=conj, t_um=t, corr=float(corr)))
    return sorted(out, key=lambda r: -r['corr'])


def amplitude_corr(T_real, T_sim, dx, radius_um):
    from scipy.ndimage import gaussian_filter
    n = T_real.shape[0]
    m = np.hypot(*(np.mgrid[:n, :n] - n // 2) * dx) < radius_um
    a = gaussian_filter(np.abs(T_real), 1.0) - 1
    return float(np.corrcoef(a[m], (np.abs(T_sim) - 1)[m])[0, 1])


def carrier_bins(geom, shape):
    ny, nx = shape
    bx = int(geom['cols'][geom['patch_center_col'] - 1])
    by = int(geom['rows'][geom['patch_center_row'] - 1])
    fx = bx if bx < nx / 2 else bx - nx
    fy = by if by < ny / 2 else by - ny
    return -float(fx), -float(fy)


def reference_field(shape, carrier, amplitude=1.0):
    ny, nx = shape
    m, n = np.mgrid[:ny, :nx]
    return amplitude * np.exp(2j * np.pi * (carrier[0] * n / nx + carrier[1] * m / ny))
