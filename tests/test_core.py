# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from fdtd_meep.acquisition import fourier_geometry
from fdtd_meep.camera import DIHEDRAL, carrier_bins, reference_field, to_camera_grid
from fdtd_meep.config import RunConfig, SimulationParams
from fdtd_meep.illuminations import generate_illuminations, match_to_holograms, read_illumination_list
from fdtd_meep.mie import angular_spectrum
from fdtd_meep.objects import (auto_crop, medium_index_from_file, read_voxel_h5, registered_gt_crop,
                               rotate_point, rotate_volume_zyx, sphere_volume, write_voxel_h5)
from fdtd_meep.scattered import cell_layout, discrete_plane_wave
from fdtd_meep.specular import analyze_hologram, fps_order, window2d

ROOT = Path(__file__).resolve().parents[1]
CFG = {'NA': 1.3, 'N0': 1.5, 'LAMBDA': 632.8e-9, 'F_TUBE': 0.15, 'F_OBJ': 0.0018,
       'TPCAM': 8e-6, 'RF': 0.75, 'DIM_ROI': 1024, 'CIRCLE_CX': 192, 'CIRCLE_CY': 192}


class MieTests(unittest.TestCase):
    def test_weak_sphere_matches_born(self):
        lam, nb, ns, R = 0.6328, 1.5, 1.503, 0.3
        k = 2 * np.pi * nb / lam
        A = angular_spectrum(np.array([0.0]), np.array([0.0]), [0, 0, k], [1, 0, 0], ns, nb, R, lam)[0, 0]
        f = A * k / (2j * np.pi)
        f_born = k ** 2 / (4 * np.pi) * ((ns / nb) ** 2 - 1) * 4 / 3 * np.pi * R ** 3
        self.assertAlmostEqual(abs(f / f_born), 1.0, delta=0.01)
        self.assertLess(abs(np.angle(f / f_born)), 0.05)


class ScatteredTests(unittest.TestCase):
    def test_discrete_plane_wave(self):
        k = 2 * np.pi * 1.5 / 0.633
        kz = k * np.sqrt(1 - 0.4 ** 2 - 0.2 ** 2)
        errs = []
        for res in (14, 28, 56):
            pw = discrete_plane_wave(res, 0.633, 1.5, 0.4, -0.2)
            self.assertLess(abs(np.dot(pw['pol'], pw['ktilde'])), 1e-12)
            errs.append(abs(pw['kappa'] - kz))
        self.assertLess(errs[1], errs[0] / 3)
        self.assertLess(errs[2], errs[1] / 3)

    def test_cell_layout_even_and_extended_toward_tilt(self):
        p = SimulationParams(resolution=12)
        for kx, ky in ((0.0, 0.0), (0.5, -0.3)):
            pw = discrete_plane_wave(12, 0.6328, 1.5, kx, ky)
            lay = cell_layout((40, 50, 60), [0.1] * 3, (0.0, 0.0, 0.0), 2.5, pw, p)
            cells = np.round(lay['cell'] * 12).astype(int)
            self.assertTrue(np.all(cells % 2 == 0))
            if kx:
                self.assertGreater(lay['center'][0], 0)
                self.assertLess(lay['center'][1], 0)


    def test_cell_layout_reflection(self):
        p = SimulationParams(resolution=12, pml_far_um=2.0)
        pw = discrete_plane_wave(12, 0.6328, 1.5, 0.5, -0.3)
        t = cell_layout((40, 50, 60), [0.1] * 3, (0.0, 0.0, 0.0), 2.5, pw, p)
        r = cell_layout((40, 50, 60), [0.1] * 3, (0.0, 0.0, 0.0), -2.5, pw, p)
        self.assertTrue(r['reflection'] and not t['reflection'])
        np.testing.assert_allclose(np.sign(r['walk']), np.sign(t['walk']))
        self.assertLess(r['center'][2] - r['cell'][2] / 2, -2.5 - 1.0)
        self.assertEqual(r['pml'][1][2], 2.0)
        self.assertTrue(np.all(np.round(r['cell'] * 12).astype(int) % 2 == 0))
        with self.assertRaises(ValueError):
            cell_layout((40, 50, 60), [0.1] * 3, (0.0, 0.0, 0.0), 0.5, pw, p)


class ObjectTests(unittest.TestCase):
    def test_voxel_roundtrip_and_rotation(self):
        n = 1.5 + np.random.default_rng(0).random((4, 5, 6)) * 0.01
        with tempfile.TemporaryDirectory() as tmp:
            write_voxel_h5(Path(tmp) / 'o.h5', n, 0.1)
            m, step = read_voxel_h5(Path(tmp) / 'o.h5')
        np.testing.assert_allclose(m, n, rtol=1e-6)
        np.testing.assert_allclose(step, 0.1, rtol=1e-12)
        np.testing.assert_array_equal(rotate_volume_zyx(n, 'y'), n[::-1, :, ::-1])
        self.assertEqual(rotate_point((1, 2, 3), 'y'), (-1, 2, -3))

    def test_sphere_volume(self):
        n = sphere_volume(1.0, 0.03, 1.5, 0.1)
        self.assertAlmostEqual(float(n.max()), 1.53, places=3)
        self.assertAlmostEqual(float(n[0, 0, 0]), 1.5, places=6)


class CameraTests(unittest.TestCase):
    def test_mapping(self):
        x = y = np.arange(-40, 41) * 0.1
        X, Y = np.meshgrid(x, y)
        S = np.exp(-((X - 1.0) ** 2 + Y ** 2) / 0.1)
        for name, exp in (('id', (1.0, 0.0)), ('rot90', (0.0, 1.0)), ('flip_x', (-1.0, 0.0))):
            T = to_camera_grid(S, x, y, DIHEDRAL[name], (0.5, -0.5), 101, 0.1) - 1
            iy, ix = np.unravel_index(np.argmax(np.abs(T)), T.shape)
            self.assertAlmostEqual((ix - 50) * 0.1, exp[0] + 0.5, places=6)
            self.assertAlmostEqual((iy - 50) * 0.1, exp[1] - 0.5, places=6)

    def test_hologram_peak_round_trip(self):
        shape = (1024, 1024)
        geom = fourier_geometry(CFG, shape)
        self.assertEqual(carrier_bins(geom, shape), (320.0, 320.0))
        dk = geom['dk'] * 1e-6
        c = (np.arange(1024) - 512) * geom['dx_obj'] * 1e6
        X, Y = np.meshgrid(c, c)
        R = reference_field(shape, carrier_bins(geom, shape), 1.1)
        for drow, dcol in ((0, 0), (0, 40), (60, -60)):
            I = np.abs(np.exp(1j * (dcol * dk * X + drow * dk * Y)) + R) ** 2
            r = analyze_hologram(I, geom, window2d(1024, 0.1))
            self.assertAlmostEqual(r['drow_subpix'], drow, places=6)
            self.assertAlmostEqual(r['dcol_subpix'], dcol, places=6)


class SelectionTests(unittest.TestCase):
    def test_fps_starts_at_center_and_covers(self):
        rng = np.random.default_rng(1)
        k = rng.uniform(-0.7, 0.7, (200, 2))
        k[17] = (0.001, -0.002)
        order = fps_order(k[:, 0], k[:, 1])
        self.assertEqual(order[0], 17)
        self.assertEqual(sorted(order), list(range(200)))
        first = k[order[:16]]
        self.assertGreater(np.ptp(first[:, 0]), 1.0)


class IlluminationListTests(unittest.TestCase):
    def write(self, text):
        tmp = tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False)
        tmp.write(text)
        tmp.close()
        self.addCleanup(Path(tmp.name).unlink)
        return tmp.name

    def test_kx_ky_and_angles_agree(self):
        a = read_illumination_list(self.write('# c\nname,kx,ky\np0,0.3,0.4\n'))
        b = read_illumination_list(self.write('theta_deg,phi_deg\n30,53.13010235\n'))
        self.assertEqual(a[0]['name'], 'p0')
        self.assertEqual(b[0]['name'], 'i000')
        self.assertAlmostEqual(a[0]['theta_deg'], 30.0, places=6)
        self.assertAlmostEqual(b[0]['kx'], 0.3, places=6)
        self.assertAlmostEqual(b[0]['ky'], 0.4, places=6)

    def test_selection_csv_is_accepted(self):
        rows = read_illumination_list(self.write(
            'index,filename,kx_norm_raw,ky_norm_raw,new_filename\n5,i005.pgm,0.1,0.2,i000.pgm\n'))
        self.assertEqual((rows[0]['name'], rows[0]['filename'], rows[0]['kx']), ('i000', 'i005.pgm', 0.1))

    def test_non_propagating_is_rejected(self):
        with self.assertRaises(ValueError):
            read_illumination_list(self.write('kx,ky\n0.8,0.7\n'))

    def test_match_by_name_and_nearest_k(self):
        rows = [dict(filename=f'i{j:03d}.pgm', kx_norm_raw=0.1 * j, ky_norm_raw=0.0) for j in range(5)]
        picks = match_to_holograms(rows, [dict(filename='i003.pgm'), dict(kx=0.29, ky=0.0),
                                          dict(kx=0.12, ky=0.01)])
        self.assertEqual(picks, [3, 2, 1])

    def test_generated_lists(self):
        for pattern in ('spiral', 'rings'):
            rows = generate_illuminations(20, 40.0, pattern)
            self.assertEqual(len(rows), 20)
            self.assertAlmostEqual(rows[0]['theta_deg'], 0.0)
            self.assertLessEqual(max(r['theta_deg'] for r in rows), 40.0 + 1e-9)


class RegisteredGTTests(unittest.TestCase):
    def make_gt(self, tmp, n_medium_key='n_medium', value=1.5):
        import h5py
        dn = np.zeros((40, 30, 20), np.float32)          # vista h5py de V(y,x,z): [z, x, y]
        dn[10:20, 8:18, 5:12] = 0.03
        path = Path(tmp) / 'gt.mat'
        with h5py.File(path, 'w') as f:
            f['phantom_gt_dn'] = dn
            f['meta/rec_voxel_um_xyz'] = np.array([0.1, 0.1, 0.1])
            f[f'meta/{n_medium_key}'] = np.array([value])
        return path

    def test_medium_from_file_without_correction(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.make_gt(tmp, 'phantom_n_immersion', 1.512)
            self.assertAlmostEqual(medium_index_from_file(path), 1.512)
            n, _, meta = registered_gt_crop(path, [[0, 30], [0, 20], [0, 40]])
            self.assertAlmostEqual(float(n.min()), 1.512, places=6)
            self.assertAlmostEqual(float(n.max()), 1.542, places=6)
            n2, _, _ = registered_gt_crop(path, [[0, 30], [0, 20], [0, 40]], n_medium=1.5)
            self.assertAlmostEqual(float(n2.max()), 1.53, places=6)

    def test_auto_crop_is_support_plus_margin(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = self.make_gt(tmp)
            self.assertEqual(auto_crop(path, pad_um=0.2), [[6, 20], [3, 14], [8, 22]])
            n, _, meta = registered_gt_crop(path, None, pad_um=0.2)
            self.assertEqual(n.shape, (14, 11, 14))
            self.assertEqual(meta['border_max_dn'], 0.0)


class ConfigTests(unittest.TestCase):
    def test_example_config_loads(self):
        cfg = RunConfig.load(ROOT / 'inputs' / 'config' / 'phantom_4pi.json')
        self.assertEqual(set(cfg.views), {'bottom', 'top'})
        self.assertEqual(cfg.simulation.resolution, 12)

    def test_free_config_loads(self):
        cfg = RunConfig.load(ROOT / 'inputs' / 'config' / 'sphere_free.json')
        self.assertEqual(cfg.views, {})
        self.assertEqual(cfg.object.source, 'sphere')
        self.assertEqual(len(read_illumination_list(cfg.illuminations_file())), 5)

    def test_unknown_parameter_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / 'c.json'
            p.write_text(json.dumps({'views': {}, 'simulation': {'resolucion': 12}}))
            with self.assertRaises(ValueError):
                RunConfig.load(p)


if __name__ == '__main__':
    unittest.main()
