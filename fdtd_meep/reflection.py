# FDTD_meep v1.1.0 | 2026-08-27| Carlos Chacón

import csv
import json
import math
import re
import shutil
import time

import numpy as np
from PIL import Image

from .acquisition import list_holograms, read_pgm
from .camera import DIHEDRAL, carrier_bins, reference_field, to_camera_grid
from .optics import focus_field
from .specular import analyze_hologram, read_selection, window2d


def variant_dir(vdir, r_spec):
    return vdir / ('sim_subset' if r_spec == 0 else f'sim_subset_spec{r_spec:g}')


def camera_field(pipe, npz, cal_c, acq):
    """Campo retrodispersado en la cámara opuesta, relativo a la amplitud incidente.

    La simulación se hace en el marco de la iluminación con el sensor bajo el objeto. La cámara
    que registra la reflexión es la de la otra vista, calibrada con el objeto girado 180°
    alrededor de y: x -> -x, z -> -z. En ese marco la retrodispersión viaja hacia +z, como el
    campo transmitido de esa vista, y se aplica su calibración (M, t, foco) sin cambios.
    """
    S, x, y, d = pipe.relative_field(npz)
    xc, Sc = -x[::-1], S[:, ::-1]
    kt = np.array([-float(d['kx']), float(d['ky'])])
    zc = -float(d['z'])
    Sf = focus_field(Sc, xc, y, kt[0], kt[1], cal_c['z_focus_um'] - zc, acq.NA,
                     pipe.wavelength, pipe.n_medium)
    M = np.array(DIHEDRAL[cal_c['M']], float)
    T = to_camera_grid(Sf, xc, y, M, cal_c['t_um'], acq.shape[0], acq.dx_cam_um, offset=0.0)
    return T, M @ kt


def run_reflection(pipe, rview, limit=None):
    rc = pipe.cfg.reflection
    spec = rc.views[rview]
    calib = pipe.calibration()
    cal_i, cal_c = calib[spec['illumination']]['best'], calib[spec['camera']]['best']
    if cal_i['orientation'] == cal_c['orientation']:
        raise ValueError(f'{rview}: la iluminación ({spec["illumination"]}) y la cámara '
                         f'({spec["camera"]}) deben entrar por caras opuestas del objeto.')
    acq = pipe.acq[spec['camera']]
    rows = read_selection(pipe.view_dir(spec['illumination']))
    rows = rows[:limit] if limit else rows
    vdir = pipe.view_dir(rview)
    fields_dir = vdir / 'fields'
    fields_dir.mkdir(parents=True, exist_ok=True)

    walls = {}
    for j, row in enumerate(rows):
        name = row['new_filename']
        k_obj = pipe.object_k(cal_i, row)
        npz = fields_dir / name.replace('.pgm', '.npz')
        if npz.exists():
            continue
        lay = pipe.layout(cal_i['orientation'], float(k_obj[0]), float(k_obj[1]), reflection=True)
        pipe.check_memory(lay)
        print(f"[{rview} {j + 1}/{len(rows)}] {name} (iluminación {spec['illumination']} "
              f"{row['filename']}) θ={float(row['theta_deg']):.1f}° "
              f"celdas={lay['n_cells'] / 1e6:.1f} M (~{lay['memory_gb']:.1f} GB)", flush=True)
        t0 = time.time()
        pipe.run_sim(cal_i['orientation'], float(k_obj[0]), float(k_obj[1]), npz, reflection=True)
        walls[name] = time.time() - t0

    k0 = 2 * math.pi * pipe.n_medium / pipe.wavelength
    sgn = -1.0 if cal_c['conj'] else 1.0
    cams = []
    for row in rows:
        npz = fields_dir / row['new_filename'].replace('.pgm', '.npz')
        T, kc = camera_field(pipe, npz, cal_c, acq)
        cams.append((row, T.astype(np.complex64), kc))

    # Un solo brazo de referencia y una sola exposición por vista, como en el montaje.
    level = float(np.median([np.percentile(np.abs(T), 99.9) for _, T, _ in cams]))
    a_ref = rc.reference_ratio * max([level] + list(rc.parasitic_specular))
    target = read_pgm(list_holograms(acq.folder)[0]).mean()
    scale = target / (a_ref ** 2 + float(np.median([np.mean(np.abs(T) ** 2) for _, T, _ in cams])))
    fx, fy = carrier_bins(acq.geom, acq.shape)
    R = reference_field(acq.shape, (-fx, -fy) if cal_c['conj'] else (fx, fy), a_ref)
    c = (np.arange(acq.shape[0]) - acq.shape[0] // 2) * acq.dx_cam_um
    Xc, Yc = np.meshgrid(c, c)
    window = window2d(acq.shape[1], pipe.cfg.selection.tukey)
    geom = acq.geom
    t = cal_c['t_um']

    variants = [0.0] + [float(r) for r in rc.parasitic_specular if r]
    records = {r: [] for r in variants}
    cam_dir = vdir / 'camera_fields'
    cam_dir.mkdir(exist_ok=True)
    for r_spec in variants:
        variant_dir(vdir, r_spec).mkdir(exist_ok=True)
    for row, T, kc in cams:
        name = row['new_filename']
        k_raw = sgn * kc / k0
        drow_s, dcol_s = k_raw[1] * geom['REwald'], k_raw[0] * geom['REwald']
        ramp = np.exp(1j * (kc[0] * (Xc - t[0]) + kc[1] * (Yc - t[1])))
        U = T * ramp
        cross = 2 * np.real(U * np.conj(R))
        P = np.fft.fft2(window * cross)[np.ix_(geom['rows'], geom['cols'])] / (a_ref * cross.size)
        np.savez_compressed(cam_dir / name.replace('.pgm', '.npz'), spectrum_patch=P.astype(np.complex64),
                            drow_specular=drow_s, dcol_specular=dcol_s,
                            kx_norm_raw_specular=k_raw[0], ky_norm_raw_specular=k_raw[1],
                            patch_center=geom['patch_center_row'], fmaxHolo=geom['fmaxHolo'],
                            REwald=geom['REwald'], dx_patch_um=geom['FOV'] * 1e6 / P.shape[0],
                            a_ref=a_ref, gray_scale=scale)
        npz = np.load(fields_dir / name.replace('.pgm', '.npz'))
        for r_spec in variants:
            I = np.abs(U + r_spec * ramp + R) ** 2 * scale
            Iq = np.clip(np.round(I), 0, 255)
            Image.fromarray(Iq.astype(np.uint8), mode='L').save(variant_dir(vdir, r_spec) / name)
            pk = analyze_hologram(Iq, geom, window)
            records[r_spec].append(dict(
                new_filename=name, illumination_view=spec['illumination'], original=row['filename'],
                theta_inc_deg=float(row['theta_deg']),
                kx_norm_raw_illum=float(row['kx_norm_raw']), ky_norm_raw_illum=float(row['ky_norm_raw']),
                kx_norm_raw_specular=float(k_raw[0]), ky_norm_raw_specular=float(k_raw[1]),
                drow_specular=float(drow_s), dcol_specular=float(dcol_s),
                drow_peak=pk['drow_subpix'], dcol_peak=pk['dcol_subpix'],
                peak_to_specular_bins=math.hypot(pk['drow_subpix'] - drow_s, pk['dcol_subpix'] - dcol_s),
                peak_snr=pk['peak_snr'], parasitic_specular=r_spec,
                back_rms=float(np.sqrt(np.mean(np.abs(T) ** 2))), back_max=float(np.abs(T).max()),
                saturated_fraction=float((I > 255).mean()),
                meep_time=float(npz['meep_time']), n_cells_M=float(npz['n_cells']) / 1e6,
                wall_s=walls.get(name, float('nan'))))

    sel_txt = pipe.view_dir(spec['illumination']) / 'selection' / 'seleccion.txt'
    for r_spec in variants:
        out = variant_dir(vdir, r_spec)
        ref = np.clip(np.round(np.full(acq.shape, a_ref ** 2 * scale)), 0, 255)
        Image.fromarray(ref.astype(np.uint8), mode='L').save(out / 'Intensite_ref.pgm')
        txt = (acq.folder / 'config_manip.txt').read_text()
        (out / 'config_manip.txt').write_text(re.sub(r'(?m)^NB_HOLO\s+\d+', f'NB_HOLO {len(rows)}', txt))
        if (acq.folder / 'recon.txt').exists():
            shutil.copy(acq.folder / 'recon.txt', out / 'recon.txt')
        if sel_txt.exists():
            shutil.copy(sel_txt, out / 'seleccion.txt')
        with open(out / 'sim_report.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(records[r_spec][0]))
            w.writeheader()
            w.writerows(records[r_spec])
        (out / 'sim_parameters.json').write_text(json.dumps(dict(
            reflection_view=rview, illumination_view=spec['illumination'], camera_view=spec['camera'],
            calibration_illumination=cal_i, calibration_camera=cal_c, a_ref=a_ref,
            back_level_p999=level, reference_ratio=rc.reference_ratio, parasitic_specular=r_spec,
            gray_scale=scale, simulation=vars(pipe.cfg.simulation), n_medium=pipe.n_medium,
            wavelength_um=pipe.wavelength), indent=2, default=str))
    figure(vdir / f'{rview}_resumen.png', rview, spec, cams, records[0.0], vdir, acq)
    levels_figure(vdir / f'{rview}_niveles.png', rview, records)
    return records


def levels_figure(path, rview, records):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    base = records[0.0]
    th = [r['theta_inc_deg'] for r in base]
    ax[0].semilogy(th, [r['back_rms'] for r in base], 'o', label='RMS |U|')
    ax[0].semilogy(th, [r['back_max'] for r in base], 's', mfc='none', label='máx |U|')
    for r_spec in records:
        if r_spec:
            ax[0].axhline(r_spec, color='C3', ls=':', label=f'especular parásito r = {r_spec:g}')
    ax[0].set_xlabel('θ de iluminación (°)')
    ax[0].set_ylabel('amplitud / |E_inc|')
    ax[0].set_title(f'{rview}: nivel de la retrodispersión en la cámara', fontsize=9)
    ax[0].legend(fontsize=7)
    for i, (r_spec, recs) in enumerate(records.items()):
        ax[1].plot(th, [max(r['peak_to_specular_bins'], 1e-2) for r in recs], 'o^'[i % 2],
                   label='sin especular' if r_spec == 0 else f'r = {r_spec:g}')
    ax[1].set_yscale('log')
    ax[1].axhline(1, color='0.5', ls=':')
    ax[1].set_xlabel('θ de iluminación (°)')
    ax[1].set_ylabel('distancia máximo del patch - especular (bins)')
    ax[1].set_title('¿El máximo del espectro marca el especular?', fontsize=9)
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def figure(path, rview, spec, cams, records, vdir, acq, targets_deg=(0.0, 12.0, 25.0, 40.0)):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    geom = acq.geom
    picks = []
    for target in targets_deg:
        j = min(range(len(records)), key=lambda i: abs(records[i]['theta_inc_deg'] - target))
        if j not in picks:
            picks.append(j)
    fig, ax = plt.subplots(4, len(picks), figsize=(4.2 * len(picks), 16.5), squeeze=False)
    n = acq.shape[0]
    h = 128
    cc = slice(n // 2 - h, n // 2 + h)
    e = [-h * acq.dx_cam_um, h * acq.dx_cam_um] * 2
    pc = geom['patch_center_row'] - 1
    for col, j in enumerate(picks):
        rec, (row, T, kc) = records[j], cams[j]
        name = rec['new_filename']
        holo = read_pgm(variant_dir(vdir, 0.0) / name)
        P = np.load(vdir / 'camera_fields' / name.replace('.pgm', '.npz'))['spectrum_patch']
        title = f"{rview} {name} (ilum. {spec['illumination']} {rec['original']}, θ={rec['theta_inc_deg']:.0f}°)"
        ax[0, col].imshow(holo[cc, cc], cmap='gray', origin='lower', extent=e)
        ax[0, col].set_title(f'{title}\nholograma (8 bits)', fontsize=7)
        L = np.log10(np.abs(P) + 1e-12)
        ax[1, col].imshow(L, cmap='magma', origin='lower', vmin=L.max() - 4, vmax=L.max())
        ax[1, col].add_patch(plt.Circle((pc, pc), geom['fmaxHolo'], fill=False, ec='w', ls=':', lw=0.8))
        ax[1, col].plot(pc + rec['dcol_specular'], pc + rec['drow_specular'], 'c+', ms=14, mew=1.5,
                        label='especular esperado')
        ax[1, col].plot(pc + rec['dcol_peak'], pc + rec['drow_peak'], 'gx', ms=10, mew=1.5,
                        label='máximo del patch')
        ax[1, col].legend(fontsize=6, loc='lower right')
        ax[1, col].set_title(f"log10|espectro| orden +1 (4 décadas)\nmáximo a {rec['peak_to_specular_bins']:.1f} bins "
                             f"del especular", fontsize=7)
        A = np.abs(T[cc, cc])
        im = ax[2, col].imshow(A, cmap='gray', origin='lower', extent=e, vmin=0, vmax=np.percentile(A, 99.9))
        fig.colorbar(im, ax=ax[2, col], fraction=0.046)
        ax[2, col].set_title('|U| retrodispersado / |E_inc|', fontsize=7)
        ax[3, col].imshow(np.angle(T[cc, cc]), cmap='twilight', origin='lower', extent=e, vmin=-np.pi, vmax=np.pi)
        ax[3, col].set_title('fase de U (sin portadora)', fontsize=7)
        for a in ax[:, col]:
            a.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)
