# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón
# 4Pi pipeline. v1.2: runs without views (free mode), sphere object, selection from a list.

import csv
import json
import math
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image

from .acquisition import Acquisition, read_pgm
from .camera import (DIHEDRAL, amplitude_corr, carrier_bins, demodulate, match,
                     reference_field, remove_background, to_camera_grid)
from .illuminations import read_illumination_list
from .objects import (medium_index_from_file, read_voxel_h5, registered_gt_crop, rotate_point,
                      rotate_volume_zyx, sphere_volume, write_voxel_h5)
from .optics import focus_field
from .specular import analyze_hologram, read_selection, select_view, window2d

ORIENTATIONS = {'zmin': None, 'zmax': 'y'}


class Pipeline:
    def __init__(self, config, config_path):
        self.cfg = config
        self.config_path = Path(config_path).resolve()
        self.out = config.output
        self.acq = {v: Acquisition(p) for v, p in config.views.items()}
        first = next(iter(self.acq.values()), None)
        # Índice del medio: el de la configuración, si no el que indica el objeto, si no N0
        from_file = (medium_index_from_file(config.object_file())
                     if config.object.source in ('registered_gt', 'voxel_h5') else None)
        self.n_medium = config.object.n_medium or from_file or (first.n0 if first else None)
        self.n_medium_source = ('config' if config.object.n_medium else 'objeto' if from_file
                                else 'config_manip N0')
        self.wavelength = first.wavelength_um if first else config.wavelength_um
        if self.n_medium is None or self.wavelength is None:
            raise ValueError('Sin vistas, object.n_medium y wavelength_um son obligatorios.')

    def view_dir(self, view):
        return self.out / view

    def object_path(self):
        return self.out / 'object' / 'object.h5'

    def prepare(self):
        o = self.cfg.object
        path = self.object_path()
        print(f'n_medio = {self.n_medium:g} ({self.n_medium_source})', flush=True)
        for view, acq in self.acq.items():
            if abs(acq.n0 - self.n_medium) > 1e-4:
                print(f'  AVISO: N0 de config_manip ({view}) = {acq.n0:g} es distinto de n_medio.',
                      flush=True)
        if o.source == 'registered_gt':
            n, voxel, meta = registered_gt_crop(self.cfg.object_file(), o.crop_xyz, self.n_medium,
                                                o.n_immersion_gt, o.soft_mask_dn, o.crop_pad_um)
            print(f"  recorte {meta['crop_xyz']} ({n.shape[2]}x{n.shape[1]}x{n.shape[0]} vóxeles), "
                  f"Δn máx. en el borde = {meta['border_max_dn']:.2g}", flush=True)
            if meta['border_max_dn'] > 1e-3:
                print('  AVISO: el objeto toca el borde del recorte; amplía crop_xyz.', flush=True)
            write_voxel_h5(path, n, voxel, meta)
        elif o.source == 'voxel_h5':
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(self.cfg.object_file(), path)
        elif o.source == 'sphere':
            n = sphere_volume(o.radius_um, o.dn, self.n_medium, o.voxel_um, edge_um=o.edge_um)
            write_voxel_h5(path, n, o.voxel_um, dict(radius_um=o.radius_um, dn=o.dn,
                                                     n_medium=self.n_medium, edge_um=o.edge_um))
        else:
            raise ValueError(f'object.source desconocido: {o.source}')
        return path

    def load_object(self, orientation, reflection=False):
        n, step = read_voxel_h5(self.object_path())
        axis = ORIENTATIONS[orientation]
        n = rotate_volume_zyx(n, axis)
        center = rotate_point(self.cfg.object.center_um, axis)
        offset = n.shape[0] * step[2] / 2 + self.cfg.simulation.sensor_gap_um
        sensor_z = center[2] - offset if reflection else center[2] + offset
        return np.ascontiguousarray(n), step, center, sensor_z

    def select(self):
        out = {}
        for view, folder in self.cfg.views.items():
            lst = self.cfg.selection_list(view)
            chosen = select_view(view, folder, self.view_dir(view), self.cfg.selection.n,
                                 self.cfg.selection.tukey,
                                 read_illumination_list(lst) if lst else None)
            out[view] = chosen
        return out

    def simulate_one(self, orientation, kx, ky, out_npz, reflection=False, slices=(), z_pad_um=0.0,
                     slice_center=None):
        import dataclasses
        import meep as mp
        from .scattered import simulate_scattered
        n, step, center, sensor_z = self.load_object(orientation, reflection)
        params = self.cfg.simulation
        if z_pad_um:
            params = dataclasses.replace(params, z_margin_um=params.z_margin_um + z_pad_um)
        r = simulate_scattered(n, step, center, self.n_medium, self.wavelength, kx, ky, sensor_z,
                               params, slices=slices, slice_center=slice_center)
        extra = {}
        for plane, sl in r['slices'].items():
            extra.update({f'slice_{plane}_{a}': sl[a] for a in 'xyz'})
            extra.update({f'slice_{plane}_Es_{c}': v for c, v in sl['Es'].items()})
            extra.update({f'slice_{plane}_Einc_{c}': v for c, v in sl['Einc'].items()})
        if mp.am_master():
            pw, lay = r['pw'], r['layout']
            np.savez_compressed(out_npz, **{f'Es_{c}': v for c, v in r['fields'].items()},
                                **{f'Einc_{c}': v for c, v in r['incident'].items()},
                                x=r['x'], y=r['y'], z=r['z'], kx=pw['kx'], ky=pw['ky'],
                                kappa=pw['kappa'], pol=pw['pol'], meep_time=r['meep_time'],
                                cell=lay['cell'], cell_center=lay['center'], walkoff=lay['walk'],
                                n_cells=lay['n_cells'], kx_norm=kx, ky_norm=ky,
                                orientation=orientation, reflection=reflection,
                                resolution=self.cfg.simulation.resolution,
                                slice_planes=np.array(list(r['slices']), dtype='U2'),
                                z_pad_um=z_pad_um, **extra)

    def run_sim(self, orientation, kx, ky, out_npz, reflection=False, extra_args=()):
        out_npz = Path(out_npz)
        out_npz.parent.mkdir(parents=True, exist_ok=True)
        nproc = self.cfg.simulation.mpi_processes
        cmd = [sys.executable, '-m', 'fdtd_meep', 'sim', '--config', str(self.config_path),
               '--orientation', orientation, f'--k={kx},{ky}', '--out', str(out_npz)]
        if reflection:
            cmd.append('--reflection')
        cmd += list(extra_args)
        if nproc > 1:
            cmd = [str(Path(sys.executable).parent / 'mpirun'), '-np', str(nproc)] + cmd
        with open(out_npz.with_suffix('.log'), 'w') as f:
            subprocess.run(cmd, cwd=Path(__file__).resolve().parents[1], stdout=f,
                           stderr=subprocess.STDOUT, check=True)

    def layout(self, orientation, kx, ky, reflection=False):
        from .scattered import cell_layout, discrete_plane_wave
        n, step, center, sensor_z = self.load_object(orientation, reflection)
        p = self.cfg.simulation
        pw = discrete_plane_wave(p.resolution, self.wavelength, self.n_medium, kx, ky, p.polarization)
        return cell_layout(n.shape, step, center, sensor_z, pw, p)

    def check_memory(self, lay):
        limit = self.cfg.simulation.max_memory_gb
        if limit and lay['memory_gb'] > limit:
            raise MemoryError(f"La celda estimada ({lay['n_cells'] / 1e6:.1f} M celdas, "
                              f"{lay['memory_gb']:.1f} GB) supera max_memory_gb = {limit}.")

    @staticmethod
    def relative_field(npz, component='Ex'):
        d = np.load(npz)
        return d[f'Es_{component}'] / d[f'Einc_{component}'], d['x'], d['y'], d

    def real_center_field(self, view):
        acq = self.acq[view]
        s = read_selection(self.view_dir(view))[0]
        U, dx = demodulate(read_pgm(acq.folder / s['filename']), acq.geom,
                           float(s['drow_subpix']), float(s['dcol_subpix']), self.cfg.selection.tukey)
        c = self.cfg.calibration
        return remove_background(U, dx, c.background_annulus_um, c.background_poly_order), dx, s

    def calibrate(self):
        c = self.cfg.calibration
        cal_dir = self.out / 'calibration'
        cal_dir.mkdir(parents=True, exist_ok=True)
        for orientation in ORIENTATIONS:
            npz = cal_dir / f'center_{orientation}.npz'
            if not npz.exists():
                self.check_memory(self.layout(orientation, 0.0, 0.0))
                print(f'Simulando iluminación normal, entrada por la cara {orientation}...', flush=True)
                self.run_sim(orientation, 0.0, 0.0, npz)
        report = {}
        for view, acq in self.acq.items():
            T_real, dx, s = self.real_center_field(view)
            dim = T_real.shape[0]
            table = []
            for orientation in ORIENTATIONS:
                S, x, y, d = self.relative_field(cal_dir / f'center_{orientation}.npz')
                zs = float(d['z'])
                best = None
                for zf in np.arange(c.z_range_um[0], c.z_range_um[1] + 1e-9, c.z_step_um):
                    Sf = focus_field(S, x, y, 0.0, 0.0, zf - zs, acq.NA, self.wavelength, self.n_medium)
                    r = match(T_real, Sf, x, y, dim, dx)[0]
                    if best is None or r['corr'] > best['corr']:
                        best = dict(r, z_focus_um=float(zf))
                z0 = best['z_focus_um']
                for zf in np.arange(z0 - c.z_step_um, z0 + c.z_step_um + 1e-9, c.z_refine_um):
                    Sf = focus_field(S, x, y, 0.0, 0.0, zf - zs, acq.NA, self.wavelength, self.n_medium)
                    r = [q for q in match(T_real, Sf, x, y, dim, dx, [best['M']])
                         if q['conj'] == best['conj']][0]
                    if r['corr'] > best['corr']:
                        best = dict(r, z_focus_um=float(zf))
                Sf = focus_field(S, x, y, 0.0, 0.0, best['z_focus_um'] - zs, acq.NA,
                                 self.wavelength, self.n_medium)
                T = to_camera_grid(Sf, x, y, DIHEDRAL[best['M']], best['t_um'], dim, dx)
                best.update(orientation=orientation,
                            corr_amplitude=amplitude_corr(T_real, T, dx, c.match_radius_um))
                table.append(best)
                print(f"  vista {view}: entrada por {orientation}: corr={best['corr']:.3f} "
                      f"amplitud={best['corr_amplitude']:.3f} M={best['M']} conj={best['conj']} "
                      f"t={np.round(best['t_um'], 2).tolist()} foco={best['z_focus_um']:+.2f} µm",
                      flush=True)
            report[view] = {'center_hologram': s['filename'], 'all': table, 'dx_patch_um': dx}
        views = list(self.acq)
        pick = lambda v, o: [r for r in report[v]['all'] if r['orientation'] == o][0]
        if len(views) == 2:
            scores = {}
            for o in ORIENTATIONS:
                other = 'zmax' if o == 'zmin' else 'zmin'
                scores[o] = sum(pick(v, ob)['corr'] + pick(v, ob)['corr_amplitude']
                                for v, ob in ((views[0], o), (views[1], other)))
            o = max(scores, key=scores.get)
            report[views[0]]['best'] = pick(views[0], o)
            report[views[1]]['best'] = pick(views[1], 'zmax' if o == 'zmin' else 'zmin')
            report['pairing_scores'] = {f'{views[0]}<-{ob}, {views[1]}<-{"zmax" if ob == "zmin" else "zmin"}': v
                                        for ob, v in scores.items()}
        else:
            for v in views:
                report[v]['best'] = max(report[v]['all'], key=lambda r: r['corr'] + r['corr_amplitude'])
        (cal_dir / 'calibration.json').write_text(json.dumps(report, indent=2))
        self.calibration_figure(report)
        return report

    def calibration(self):
        return json.loads((self.out / 'calibration' / 'calibration.json').read_text())

    def calibration_figure(self, report):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        views = list(self.acq)
        fig, ax = plt.subplots(len(views), 4, figsize=(17, 4.3 * len(views)), squeeze=False)
        for i, view in enumerate(views):
            T_real, dx, s = self.real_center_field(view)
            b = report[view]['best']
            S, x, y, d = self.relative_field(self.out / 'calibration' / f"center_{b['orientation']}.npz")
            Sf = focus_field(S, x, y, 0.0, 0.0, b['z_focus_um'] - float(d['z']), self.acq[view].NA,
                             self.wavelength, self.n_medium)
            T = to_camera_grid(Sf, x, y, DIHEDRAL[b['M']], b['t_um'], T_real.shape[0], dx)
            T = np.conj(T) if b['conj'] else T
            h = T_real.shape[0] // 2
            c = slice(h - 90, h + 90)
            e = [-90 * dx, 90 * dx, -90 * dx, 90 * dx]
            for j, (img, title, cmap) in enumerate([
                    (np.angle(T_real[c, c]), f'{view}: fase real ({s["filename"]})', 'twilight'),
                    (np.angle(T[c, c]), f'fase simulada (entrada {b["orientation"]})', 'twilight'),
                    (np.abs(T_real[c, c]), 'amplitud real', 'gray'),
                    (np.abs(T[c, c]), 'amplitud simulada', 'gray')]):
                kw = dict(vmin=-np.pi, vmax=np.pi) if cmap == 'twilight' else dict(vmin=0.6, vmax=1.4)
                im = ax[i, j].imshow(img, origin='lower', cmap=cmap, extent=e, **kw)
                ax[i, j].set_title(title, fontsize=9)
                fig.colorbar(im, ax=ax[i, j], fraction=0.046)
        fig.tight_layout()
        fig.savefig(self.out / 'calibration' / 'calibration.png', dpi=110)
        plt.close(fig)

    def object_k(self, cal, row):
        M = np.array(DIHEDRAL[cal['M']], float)
        sgn = -1.0 if cal['conj'] else 1.0
        k_cam = sgn * np.array([float(row['kx_norm_raw']), float(row['ky_norm_raw'])])
        return M.T @ k_cam

    def estimate(self):
        cal = self.calibration() if (self.out / 'calibration' / 'calibration.json').exists() else None
        rows_out = []
        for view in self.acq:
            for row in read_selection(self.view_dir(view)):
                if cal:
                    b = cal[view]['best']
                    k, orient = self.object_k(b, row), b['orientation']
                else:
                    k, orient = (float(row['kx_norm_raw']), float(row['ky_norm_raw'])), 'zmin'
                lay = self.layout(orient, float(k[0]), float(k[1]))
                rows_out.append(dict(view=view, hologram=row['new_filename'], original=row['filename'],
                                     theta_deg=float(row['theta_deg']),
                                     cell_um=' x '.join(f'{c:.2f}' for c in lay['cell']),
                                     mcells=round(lay['n_cells'] / 1e6, 2),
                                     memory_gb=round(lay['memory_gb'], 2)))
        return rows_out

    def run_view(self, view, limit=None):
        cal = self.calibration()[view]['best']
        acq = self.acq[view]
        M = np.array(DIHEDRAL[cal['M']], float)
        rows = read_selection(self.view_dir(view))
        rows = rows[:limit] if limit else rows
        fields_dir = self.view_dir(view) / 'fields'
        out_dir = self.view_dir(view) / 'sim_subset'
        fields_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)
        ref = read_pgm(acq.folder / 'Intensite_ref.pgm').mean()
        holo_mean = read_pgm(acq.folder / rows[0]['filename']).mean()
        a_ref = math.sqrt(ref / max(holo_mean - ref, 1e-6))
        fx, fy = carrier_bins(acq.geom, acq.shape)
        carrier = (-fx, -fy) if cal['conj'] else (fx, fy)
        R = reference_field(acq.shape, carrier, a_ref)
        c = (np.arange(acq.shape[0]) - acq.shape[0] // 2) * acq.dx_cam_um
        Xc, Yc = np.meshgrid(c, c)
        window = window2d(acq.shape[1], self.cfg.selection.tukey)
        records = []
        for j, row in enumerate(rows):
            name = row['new_filename']
            k_obj = self.object_k(cal, row)
            npz = fields_dir / name.replace('.pgm', '.npz')
            t0 = time.time()
            if not npz.exists():
                lay = self.layout(cal['orientation'], float(k_obj[0]), float(k_obj[1]))
                self.check_memory(lay)
                print(f"[{view} {j + 1}/{len(rows)}] {name} ({row['filename']}) θ={float(row['theta_deg']):.1f}° "
                      f"celdas={lay['n_cells'] / 1e6:.1f} M (~{lay['memory_gb']:.1f} GB)", flush=True)
                self.run_sim(cal['orientation'], float(k_obj[0]), float(k_obj[1]), npz)
            wall = time.time() - t0
            S, x, y, d = self.relative_field(npz)
            kx, ky = float(d['kx']), float(d['ky'])
            Sf = focus_field(S, x, y, kx, ky, cal['z_focus_um'] - float(d['z']), acq.NA,
                             self.wavelength, self.n_medium)
            T = to_camera_grid(Sf, x, y, M, cal['t_um'], acq.shape[0], acq.dx_cam_um)
            kc = M @ np.array([kx, ky])
            U = T * np.exp(1j * (kc[0] * (Xc - cal['t_um'][0]) + kc[1] * (Yc - cal['t_um'][1])))
            I = np.abs(U + R) ** 2
            scale = read_pgm(acq.folder / row['filename']).mean() / I.mean()
            Iq = np.clip(np.round(I * scale), 0, 255)
            Image.fromarray(Iq.astype(np.uint8), mode='L').save(out_dir / name)
            pk = analyze_hologram(Iq, acq.geom, window)
            dp = math.hypot(pk['drow_subpix'] - float(row['drow_subpix']),
                            pk['dcol_subpix'] - float(row['dcol_subpix']))
            records.append(dict(new_filename=name, original=row['filename'],
                                theta_deg=float(row['theta_deg']),
                                kx_norm_raw=float(row['kx_norm_raw']), ky_norm_raw=float(row['ky_norm_raw']),
                                kx_norm_obj=float(k_obj[0]), ky_norm_obj=float(k_obj[1]),
                                delta_p_bins=dp, saturated_fraction=float((I * scale > 255).mean()),
                                n_cells_M=float(d['n_cells']) / 1e6 if 'n_cells' in d else float('nan'),
                                meep_time=float(d['meep_time']), wall_s=wall))
            print(f'    Δp = {dp:.3f} bins', flush=True)
            if j == 0:
                ref_img = np.clip(np.round(np.full(acq.shape, a_ref ** 2 * scale)), 0, 255)
                Image.fromarray(ref_img.astype(np.uint8), mode='L').save(out_dir / 'Intensite_ref.pgm')
        real = self.view_dir(view) / 'real_subset'
        for f in ('config_manip.txt', 'recon.txt', 'seleccion.txt'):
            if (real / f).exists():
                shutil.copy(real / f, out_dir / f)
        with open(out_dir / 'sim_report.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(records[0]))
            w.writeheader()
            w.writerows(records)
        (out_dir / 'sim_parameters.json').write_text(json.dumps(
            dict(calibration=cal, reference_ratio=a_ref ** 2, carrier_bins=list(carrier),
                 simulation=vars(self.cfg.simulation), n_medium=self.n_medium,
                 wavelength_um=self.wavelength), indent=2, default=str))
        return records

    def compare(self, targets_deg=(0.0, 12.0, 40.0)):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        views = list(self.acq)
        out_dir = self.out / 'comparison'
        out_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(4 * len(views), len(targets_deg),
                               figsize=(4.2 * len(targets_deg), 8.3 * len(views)), squeeze=False)
        rows_out = []
        for iv, view in enumerate(views):
            acq = self.acq[view]
            window = window2d(acq.shape[1], self.cfg.selection.tukey)
            sel = read_selection(self.view_dir(view))
            for j, target in enumerate(targets_deg):
                row = min(sel, key=lambda r: abs(float(r['theta_deg']) - target))
                name = row['new_filename']
                T = {}
                for kind, sub in (('real', 'real_subset'), ('sim', 'sim_subset')):
                    holo = read_pgm(self.view_dir(view) / sub / name)
                    pk = analyze_hologram(holo, acq.geom, window)
                    U, dx = demodulate(holo, acq.geom, pk['drow_subpix'], pk['dcol_subpix'])
                    T[kind] = remove_background(U, dx, self.cfg.calibration.background_annulus_um,
                                                self.cfg.calibration.background_poly_order)
                n = T['real'].shape[0]
                m = np.hypot(*(np.mgrid[:n, :n] - n // 2) * dx) < self.cfg.calibration.match_radius_um
                a, b = T['real'][m] - 1, T['sim'][m] - 1
                corr = abs(np.vdot(b, a)) / (np.linalg.norm(a) * np.linalg.norm(b))
                rows_out.append(dict(view=view, hologram=name, original=row['filename'],
                                     theta_deg=float(row['theta_deg']), corr_complex=float(corr)))
                c = slice(n // 2 - 90, n // 2 + 90)
                e = [-90 * dx, 90 * dx, -90 * dx, 90 * dx]
                for k, (img, lab, cmap, kw) in enumerate([
                        (np.angle(T['real'][c, c]), 'fase real', 'twilight', dict(vmin=-np.pi, vmax=np.pi)),
                        (np.angle(T['sim'][c, c]), 'fase simulada', 'twilight', dict(vmin=-np.pi, vmax=np.pi)),
                        (np.abs(T['real'][c, c]), 'amplitud real', 'gray', dict(vmin=0.6, vmax=1.4)),
                        (np.abs(T['sim'][c, c]), 'amplitud simulada', 'gray', dict(vmin=0.6, vmax=1.4))]):
                    a_ = ax[4 * iv + k, j]
                    a_.imshow(img, origin='lower', cmap=cmap, extent=e, **kw)
                    a_.set_title(f"{view} {name} ({row['filename']}, θ={float(row['theta_deg']):.0f}°): {lab}",
                                 fontsize=7)
                    a_.tick_params(labelsize=6)
        fig.tight_layout()
        fig.savefig(out_dir / 'comparison.png', dpi=90)
        plt.close(fig)
        (out_dir / 'comparison.json').write_text(json.dumps(rows_out, indent=2))
        return rows_out
