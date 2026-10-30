# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón

import cmath
import math

import numpy as np

COURANT = 0.5
YEE_OFFSET = {'Ex': (0.5, 0.0, 0.0), 'Ey': (0.0, 0.5, 0.0), 'Ez': (0.0, 0.0, 0.5)}
POLARIZATION = {'Ex': (1, 0, 0), 'Ey': (0, 1, 0)}
SF_CALIBRATION_BY_RES = {12: 1.190, 14: 1.056, 20: 1.070}


def sf_calibration(resolution):
    r = sorted(SF_CALIBRATION_BY_RES)
    return float(np.interp(resolution, r, [SF_CALIBRATION_BY_RES[k] for k in r]))


def discrete_plane_wave(resolution, wavelength_um, n_medium, kx_norm, ky_norm, polarization='Ex'):
    d = 1.0 / resolution
    dt = COURANT * d
    w = 2 * math.pi / wavelength_um
    k = 2 * math.pi * n_medium / wavelength_um
    kx, ky = k * kx_norm, k * ky_norm
    s2 = (d / dt) ** 2 * n_medium ** 2 * math.sin(w * dt / 2) ** 2 \
        - math.sin(kx * d / 2) ** 2 - math.sin(ky * d / 2) ** 2
    if not 0 < s2 < 1:
        raise ValueError('La iluminación no es propagante en la malla discreta.')
    kappa = 2 / d * math.asin(math.sqrt(s2))
    ktil = 2 / d * np.array([math.sin(kx * d / 2), math.sin(ky * d / 2), math.sin(kappa * d / 2)])
    p0 = np.asarray(POLARIZATION[polarization], complex)
    khat = ktil / np.linalg.norm(ktil)
    e = p0 - np.dot(khat, p0) * khat
    return dict(kx=kx, ky=ky, kappa=kappa, ktilde=ktil, pol=e / np.linalg.norm(e),
                omega=w, d=d, dt=dt)


def source_spectrum(src_time, freq, dt):
    t = np.arange(0.0, src_time.swigobj.last_time() + dt, dt)
    J = np.array([complex(src_time.swigobj.current(ti, dt)) for ti in t])
    return complex(np.sum(J * np.exp(2j * np.pi * freq * t)) * dt) / math.sqrt(2 * math.pi)


def volume_axes(shape_zyx, step, center):
    return [(np.arange(shape_zyx[2 - i]) - (shape_zyx[2 - i] - 1) / 2) * step[i] + center[i]
            for i in range(3)]


def cell_layout(shape_zyx, step, center, sensor_z, pw, params):
    d = pw['d']
    axes = volume_axes(shape_zyx, step, center)
    lo = np.array([ax[0] - st / 2 for ax, st in zip(axes, step)])
    hi = np.array([ax[-1] + st / 2 for ax, st in zip(axes, step)])
    lo = np.floor(lo / d) * d - d
    hi = np.ceil(hi / d) * d + d
    reflection = sensor_z < lo[2]
    if not reflection and sensor_z <= hi[2]:
        raise ValueError(f'sensor_z ({sensor_z}) debe estar fuera de la caja del objeto '
                         f'[{lo[2]:.3f}, {hi[2]:.3f}].')
    # Transmisión: la luz sale hacia +z; reflexión: la retrodispersión baja hacia el sensor
    # con el mismo k transversal, así que el walk-off tiene el mismo signo.
    depth = hi[2] - sensor_z if reflection else sensor_z - lo[2]
    walk = depth * np.array([pw['kx'], pw['ky']]) / pw['kappa']
    c_lo, c_hi = lo.copy(), hi.copy()
    for a in range(2):
        if walk[a] > 0:
            c_hi[a] += walk[a] + params.walkoff_margin_um
        elif walk[a] < 0:
            c_lo[a] += walk[a] - params.walkoff_margin_um
        c_lo[a] -= params.lateral_margin_um
        c_hi[a] += params.lateral_margin_um
    if reflection:
        c_lo[2] = min(lo[2] - params.z_margin_um, sensor_z - 0.5)
        c_hi[2] += params.z_margin_um
    else:
        c_lo[2] -= params.z_margin_um
        c_hi[2] = max(hi[2] + params.z_margin_um, sensor_z + 0.5)
    pml = pml_thickness(params, reflection)
    c_lo, c_hi = c_lo - pml[0], c_hi + pml[1]
    c_lo = np.floor(c_lo / d) * d
    c_hi = np.ceil(c_hi / d) * d
    c_hi = c_hi + (np.round((c_hi - c_lo) / d).astype(int) % 2 == 1) * d
    cell = c_hi - c_lo
    n_cells = float(np.prod(np.round(cell / d)))
    return dict(lo=lo, hi=hi, cell=cell, center=(c_lo + c_hi) / 2, walk=walk, pml=pml,
                reflection=bool(reflection), n_cells=n_cells,
                memory_gb=n_cells * params.bytes_per_cell / 1e9)


def pml_thickness(params, reflection):
    # (lado bajo, lado alto) por eje. En reflexión, el campo dispersado hacia delante, mucho
    # más intenso que el retrodispersado, sale por +z: allí se usa pml_far_um.
    lo = np.full(3, float(params.pml_um))
    hi = lo.copy()
    if reflection and params.pml_far_um:
        hi[2] = float(params.pml_far_um)
    return lo, hi


def delta_eps_grid(eps, axes, step, eps_b, start, shape, d):
    from scipy.interpolate import RegularGridInterpolator
    ax, ay, az = axes
    coords = [start[a] + np.arange(shape[a]) * d for a in range(3)]
    out = np.zeros(shape, dtype=np.float32)
    f = RegularGridInterpolator((az, ay, ax), eps - eps_b, bounds_error=False, fill_value=0.0)
    sel = [np.flatnonzero((coords[a] >= axis[0] - step[a] / 2) & (coords[a] <= axis[-1] + step[a] / 2))
           for a, axis in enumerate((ax, ay, az))]
    if any(len(s) == 0 for s in sel):
        return out
    X, Y = np.meshgrid(coords[0][sel[0]], coords[1][sel[1]], indexing='ij')
    for iz in sel[2]:
        zq = np.clip(coords[2][iz], az[0], az[-1])
        pts = np.stack([np.full(X.size, zq), np.clip(Y.ravel(), ay[0], ay[-1]),
                        np.clip(X.ravel(), ax[0], ax[-1])], -1)
        out[np.ix_(sel[0], sel[1], [iz])] = f(pts).reshape(X.shape + (1,)).astype(np.float32)
    return out


def _output_factor(coord, offset, k, d):
    frac = (coord / d - offset) % 1.0
    return 1.0 if min(frac, 1 - frac) < 1e-6 else math.cos(k * d / 2)


def incident_at_sensor(pw, comp, x, y, z, S):
    off = YEE_OFFSET[comp]
    ks = (pw['kx'], pw['ky'], pw['kappa'])
    f = 1.0
    for a, coord in enumerate((x[len(x) // 2], y[len(y) // 2], z)):
        f *= _output_factor(coord, off[a], ks[a], pw['d'])
    X, Y = np.meshgrid(x, y)
    return S * f * pw['pol']['xyz'.index(comp[1])] * np.exp(
        1j * (pw['kx'] * X + pw['ky'] * Y + pw['kappa'] * z))


def boundary_layers(mp, pml):
    lo, hi = pml
    if np.allclose(lo, lo[0]) and np.allclose(hi, lo[0]):
        return [mp.PML(float(lo[0]))]
    out = []
    for a, direction in enumerate((mp.X, mp.Y, mp.Z)):
        out.append(mp.PML(float(lo[a]), direction=direction, side=mp.Low))
        out.append(mp.PML(float(hi[a]), direction=direction, side=mp.High))
    return out


def incident_on_grid(pw, comp, x, y, z, S):
    # Onda incidente discreta en una malla [x, y, z] de salida DFT (cortes laterales).
    off = YEE_OFFSET[comp]
    ks = (pw['kx'], pw['ky'], pw['kappa'])
    f = [np.array([_output_factor(c, off[a], ks[a], pw['d']) for c in coord])
         for a, coord in enumerate((x, y, z))]
    return (S * pw['pol']['xyz'.index(comp[1])] * f[0][:, None, None] * f[1][None, :, None]
            * f[2][None, None, :] * np.exp(1j * (pw['kx'] * x[:, None, None] + pw['ky'] * y[None, :, None]
                                                 + pw['kappa'] * z[None, None, :])))


def simulate_scattered(n_zyx, step, center, n_medium, wavelength_um, kx_norm, ky_norm,
                       sensor_z, params, real_fields=True, slices=(), slice_center=None):
    import meep as mp
    pw = discrete_plane_wave(params.resolution, wavelength_um, n_medium, kx_norm, ky_norm,
                             params.polarization)
    d = pw['d']
    lay = cell_layout(n_zyx.shape, step, center, sensor_z, pw, params)
    lo, hi = lay['lo'], lay['hi']
    eps = np.asarray(n_zyx, float) ** 2
    eps_b = n_medium ** 2
    axes = volume_axes(n_zyx.shape, step, center)

    e_lo, e_hi = float(eps.min()), float(eps.max())
    weights = np.zeros_like(eps) if e_hi == e_lo else (eps - e_lo) / (e_hi - e_lo)
    grid = mp.MaterialGrid(mp.Vector3(*n_zyx.shape[::-1]), mp.Medium(epsilon=e_lo),
                           mp.Medium(epsilon=max(e_hi, e_lo + 1e-12)),
                           weights=np.ascontiguousarray(weights.transpose(2, 1, 0)),
                           do_averaging=False)
    size = [n_zyx.shape[2 - i] * step[i] for i in range(3)]
    geometry = [mp.Block(size=mp.Vector3(*size), center=mp.Vector3(*center), material=grid)]

    fcen = 1.0 / wavelength_um
    src_time = mp.GaussianSource(frequency=fcen, fwidth=params.pulse_fwidth * fcen)
    src_center = mp.Vector3(*((lo + hi) / 2))
    src_size = mp.Vector3(*(hi - lo))
    sources, holders = [], []
    for comp, mcomp in (('Ex', mp.Ex), ('Ey', mp.Ey), ('Ez', mp.Ez)):
        e_c = pw['pol']['xyz'.index(comp[1])]
        if abs(e_c) < 1e-9:
            continue
        off = YEE_OFFSET[comp]
        start = np.array([math.floor(lo[a] / d - off[a]) + off[a] for a in range(3)]) * d
        shape = tuple(np.round((hi - start) / d).astype(int) + 2)
        holder = {'de': delta_eps_grid(eps, axes, step, eps_b, start, shape, d)}
        holders.append(holder)
        amp0 = -1j * pw['omega'] * e_c
        kx, ky, kz = pw['kx'], pw['ky'], pw['kappa']
        cx, cy, cz = src_center.x, src_center.y, src_center.z

        def amp(v, holder=holder, start=start, amp0=amp0, cx=cx, cy=cy, cz=cz):
            de = holder['de']
            x, y, z = v.x + cx, v.y + cy, v.z + cz
            i = int(round((x - start[0]) / d))
            j = int(round((y - start[1]) / d))
            k = int(round((z - start[2]) / d))
            if i < 0 or j < 0 or k < 0 or i >= de.shape[0] or j >= de.shape[1] or k >= de.shape[2]:
                return 0j
            val = de[i, j, k]
            if val == 0.0:
                return 0j
            return amp0 * float(val) * cmath.exp(1j * (kx * x + ky * y + kz * z))
        sources.append(mp.Source(src_time, component=mcomp, center=src_center, size=src_size,
                                 amp_func=amp))

    sim = mp.Simulation(cell_size=mp.Vector3(*lay['cell']), geometry_center=mp.Vector3(*lay['center']),
                        resolution=params.resolution, geometry=geometry, sources=sources,
                        boundary_layers=boundary_layers(mp, lay['pml']),
                        default_material=mp.Medium(index=n_medium),
                        force_complex_fields=not real_fields)
    sim.init_sim()
    for holder in holders:
        holder['de'] = None
    comps, names = [mp.Ex, mp.Ey, mp.Ez], ['Ex', 'Ey', 'Ez']
    p_lo, p_hi = lay['pml']
    c_lo = lay['center'] - lay['cell'] / 2
    dft = sim.add_dft_fields(comps, fcen, 0, 1,
                             center=mp.Vector3(c_lo[0] + (p_lo[0] + lay['cell'][0] - p_hi[0]) / 2,
                                               c_lo[1] + (p_lo[1] + lay['cell'][1] - p_hi[1]) / 2,
                                               sensor_z),
                             size=mp.Vector3(lay['cell'][0] - p_lo[0] - p_hi[0],
                                             lay['cell'][1] - p_lo[1] - p_hi[1], 0))
    # Cortes laterales opcionales (XZ / YZ) de todo el dominio sin PML, para figuras de referencia.
    slice_mons = {}
    pos = (lo + hi) / 2 if slice_center is None else np.asarray(slice_center, float)
    for plane in slices:
        sc = c_lo + (p_lo + lay['cell'] - p_hi) / 2
        ss = lay['cell'] - p_lo - p_hi
        a = {'XZ': 1, 'YZ': 0}[plane]
        sc[a], ss[a] = pos[a], 0.0
        slice_mons[plane] = sim.add_dft_fields(comps, fcen, 0, 1, center=mp.Vector3(*sc),
                                               size=mp.Vector3(*ss))
    primary = params.polarization
    obj_c = (lo + hi) / 2
    probe = mp.Vector3(obj_c[0] + lay['walk'][0] / 2, obj_c[1] + lay['walk'][1] / 2, sensor_z)
    t_src_end = src_time.swigobj.last_time()
    t_max = params.max_runtime_periods / fcen
    decayed = mp.stop_when_fields_decayed(5 / fcen, comps[names.index(primary)], probe,
                                          params.field_decay_by)
    sim.run(until=[lambda s: decayed(s) and s.meep_time() > t_src_end,
                   lambda s: s.meep_time() >= t_max])

    S = source_spectrum(src_time, fcen, COURANT / params.resolution)
    if real_fields:
        S = S / 2
    fields = {name: np.asarray(sim.get_dft_array(dft, c, 0)).T for c, name in zip(comps, names)}
    xs, ys, zs, _ = sim.get_array_metadata(dft_cell=dft)
    x, y = np.asarray(xs), np.asarray(ys)
    z = float(np.asarray(zs).ravel()[0])
    cal = sf_calibration(params.resolution)
    slices_out = {}
    for plane, mon in slice_mons.items():
        sx, sy, sz, _ = sim.get_array_metadata(dft_cell=mon)
        sx, sy, sz = (np.atleast_1d(np.asarray(v, float)).ravel() for v in (sx, sy, sz))
        shape = (len(sx), len(sy), len(sz))
        slices_out[plane] = dict(
            x=sx, y=sy, z=sz,
            # Sin el factor sf_calibration: Es y Einc comparten la normalización de Meep y el
            # cociente Es/Einc ya es el físico (intensidad total ≈ 1 dentro de un objeto de fase).
            Es={name: np.asarray(sim.get_dft_array(mon, c, 0)).reshape(shape)
                for c, name in zip(comps, names)},
            Einc={name: incident_on_grid(pw, name, sx, sy, sz, S) for name in names})
    out = dict(fields={k: v * cal for k, v in fields.items()}, slices=slices_out,
               incident={name: incident_at_sensor(pw, name, x, y, z, S) for name in names},
               x=x, y=y, z=z, primary=primary, pw=pw, meep_time=sim.meep_time(),
               layout=lay, S=S)
    sim.reset_meep()
    return out
