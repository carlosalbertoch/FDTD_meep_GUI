# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Main window: parameters (left), results viewer, saved figures and console (right).

import copy
import json
from pathlib import Path

import numpy as np
from PyQt5.QtCore import Qt, QTimer, QUrl, pyqtSignal
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
                             QFileDialog, QFormLayout, QFrame, QGridLayout, QGroupBox, QHBoxLayout,
                             QLabel, QLineEdit, QListWidget, QListWidgetItem, QMainWindow,
                             QMessageBox, QPlainTextEdit, QProgressBar, QPushButton, QScrollArea,
                             QSpinBox, QSplitter, QTabWidget, QTextBrowser, QVBoxLayout, QWidget)

from fdtd_meep import __version__
from fdtd_meep.config import ROOT, RunConfig
from fdtd_meep.illuminations import (generate_illuminations, read_illumination_list,
                                     write_illumination_list)

from . import data
from .plots import Panel
from .runner import Runner

WORKFLOWS = ('4Pi acquisition (real holograms)', 'Free simulation (sphere / own object)')
SOURCES = ('sphere', 'voxel_h5', 'registered_gt')
CMAPS = ('hot', 'gray', 'viridis', 'magma', 'inferno', 'turbo')
DEFAULT_CONFIG = ROOT / 'inputs' / 'config' / 'sphere_free.json'


def dspin(value, lo, hi, step=0.1, decimals=2, special=None):
    w = QDoubleSpinBox()
    w.setRange(lo, hi)
    w.setDecimals(decimals)
    w.setSingleStep(step)
    w.setValue(value)
    if special:
        w.setSpecialValueText(special)
    return w


def ispin(value, lo, hi, special=None):
    w = QSpinBox()
    w.setRange(lo, hi)
    w.setValue(value)
    if special:
        w.setSpecialValueText(special)
    return w


def set_row_visible(form, field, visible):
    label = form.labelForField(field)
    if label:
        label.setVisible(visible)
    field.setVisible(visible)


class PathEdit(QWidget):
    changed = pyqtSignal()

    def __init__(self, kind='file', filt='', placeholder=''):
        super().__init__()
        self.kind, self.filt = kind, filt
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        self.edit.editingFinished.connect(self.changed)
        btn = QPushButton('…')
        btn.setFixedWidth(32)
        btn.setToolTip('Browse')
        btn.clicked.connect(self.browse)
        lay.addWidget(self.edit, 1)
        lay.addWidget(btn)

    def browse(self):
        start = data.resolve(self.text()) or ROOT
        start = str(start if start.exists() else ROOT)
        if self.kind == 'dir':
            p = QFileDialog.getExistingDirectory(self, 'Select folder', start)
        else:
            p, _ = QFileDialog.getOpenFileName(self, 'Select file', start, self.filt)
        if p:
            self.setText(data.relative(p))
            self.changed.emit()

    def text(self):
        return self.edit.text().strip()

    def setText(self, text):
        self.edit.setText(text or '')


class GenerateDialog(QDialog):
    """Genera una lista de iluminaciones uniforme dentro de un cono θ ≤ θmax."""

    def __init__(self, parent, default_path):
        super().__init__(parent)
        self.setWindowTitle('Generate illumination list')
        self.resize(520, 600)
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.n = ispin(20, 1, 1000)
        self.theta = dspin(35.0, 0.0, 80.0, 1.0, 1)
        self.pattern = QComboBox()
        self.pattern.addItems(['spiral', 'rings'])
        self.path = QLineEdit(default_path)
        form.addRow('Number of illuminations', self.n)
        form.addRow('Max. angle θ (deg)', self.theta)
        form.addRow('Pattern', self.pattern)
        form.addRow('Save as', self.path)
        lay.addLayout(form)
        self.preview = Panel('Preview (kx/k, ky/k)')
        lay.addWidget(self.preview, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.save)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)
        for w in (self.n, self.theta):
            w.valueChanged.connect(self.update_preview)
        self.pattern.currentIndexChanged.connect(self.update_preview)
        self.update_preview()

    def rows(self):
        return generate_illuminations(self.n.value(), self.theta.value(), self.pattern.currentText())

    def update_preview(self):
        r = self.rows()
        self.preview.scatter_sets([('generated', [q['kx'] for q in r], [q['ky'] for q in r],
                                    [q['name'] for q in r], 'sel')],
                                  rho_max=np.sin(np.radians(self.theta.value())))

    def save(self):
        path = data.resolve(self.path.text().strip())
        if path.exists() and QMessageBox.question(
                self, 'Overwrite?', f'{path.name} exists. Overwrite it?') != QMessageBox.Yes:
            return
        write_illumination_list(path, self.rows())
        self.saved_path = data.relative(path)
        self.accept()


class MainWindow(QMainWindow):
    def __init__(self, config_path=None):
        super().__init__()
        self.setWindowTitle(f'FDTD_meep 4Pi — v{__version__}')
        self.resize(1680, 1020)
        self.cfg, self.cfg_path = {}, None
        self._views_memory = {}
        self._holo_cache = {}
        self._object_key, self._object = None, None
        self._sources = {}
        self.runner = Runner(self)
        self.runner.line.connect(self.log)
        self.runner.progress.connect(self.on_progress)
        self.runner.step_done.connect(lambda _: self.refresh_viewer())
        self.runner.finished.connect(self.on_finished)
        self.preview_timer = QTimer(self, singleShot=True, interval=500)
        self.preview_timer.timeout.connect(self.update_object_panel)

        split = QSplitter(Qt.Horizontal)
        split.addWidget(self.build_left())
        split.addWidget(self.build_right())
        split.setStretchFactor(1, 1)
        split.setSizes([400, 1280])
        self.setCentralWidget(split)
        self.statusBar().showMessage('Ready.')
        self.load_config(Path(config_path) if config_path else DEFAULT_CONFIG)

    # ------------------------------------------------------------------ left panel
    def build_left(self):
        inner = QWidget()
        col = QVBoxLayout(inner)
        col.setContentsMargins(10, 10, 10, 10)
        col.setSpacing(6)
        head = QHBoxLayout()
        title = QLabel('FDTD_meep · 4Pi')
        title.setObjectName('AppTitle')
        ver = QLabel(f'v{__version__}')
        ver.setObjectName('Version')
        head.addWidget(title)
        head.addStretch()
        head.addWidget(ver)
        col.addLayout(head)
        sub = QLabel('FDTD (Meep) simulation of tomographic diffraction microscopy: '
                     'transmission and 4Pi reflection holograms.')
        sub.setObjectName('Muted')
        sub.setWordWrap(True)
        col.addWidget(sub)

        # Project
        g = QGroupBox('Project')
        f = QFormLayout(g)
        self.workflow = QComboBox()
        self.workflow.addItems(WORKFLOWS)
        self.workflow.currentIndexChanged.connect(self.on_workflow)
        f.addRow('Workflow', self.workflow)
        self.cfg_label = QLabel('—')
        self.cfg_label.setObjectName('Info')
        self.cfg_label.setWordWrap(True)
        f.addRow('Config', self.cfg_label)
        row = QHBoxLayout()
        for text, slot in (('Open…', self.open_config), ('Save', self.save_config),
                           ('Save as…', self.save_config_as)):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        f.addRow(row)
        self.name = QLineEdit()
        f.addRow('Run name', self.name)
        self.output = PathEdit('dir', placeholder='outputs/<run name>')
        self.output.changed.connect(self.refresh_viewer)
        f.addRow('Output folder', self.output)
        col.addWidget(g)

        # Real data (4Pi)
        self.g_data = QGroupBox('Input holograms (4Pi views)')
        f = QFormLayout(self.g_data)
        self.view_bottom = PathEdit('dir', placeholder='folder with config_manip.txt + iNNN.pgm')
        self.view_top = PathEdit('dir', placeholder='folder with config_manip.txt + iNNN.pgm')
        f.addRow('Bottom view', self.view_bottom)
        f.addRow('Top view', self.view_top)
        self.data_info = QLabel()
        self.data_info.setObjectName('Info')
        self.data_info.setWordWrap(True)
        f.addRow(self.data_info)
        for w in (self.view_bottom, self.view_top):
            w.changed.connect(self.update_data_info)
        col.addWidget(self.g_data)

        # Object
        g = QGroupBox('Object')
        self.f_obj = f = QFormLayout(g)
        self.obj_source = QComboBox()
        self.obj_source.addItems(SOURCES)
        self.obj_source.currentIndexChanged.connect(self.on_source)
        f.addRow('Source', self.obj_source)
        self.obj_file = PathEdit('file', 'Objects (*.h5 *.mat);;All files (*)')
        f.addRow('File', self.obj_file)
        self.crop_row = QWidget()
        crop_lay = QHBoxLayout(self.crop_row)
        crop_lay.setContentsMargins(0, 0, 0, 0)
        crop_lay.setSpacing(4)
        self.obj_crop = QLineEdit()
        self.obj_crop.setPlaceholderText('auto: object support + margin')
        self.obj_crop.setToolTip('[[x0,x1],[y0,y1],[z0,z1]]: voxel indices of the GT, 0-based, end '
                                 'excluded (MATLAB: V(y0+1:y1, x0+1:x1, z0+1:z1)). Empty = automatic.')
        crop_btn = QPushButton('Auto')
        crop_btn.setToolTip('Compute the box that contains the object (|Δn| > threshold) plus the margin')
        crop_btn.clicked.connect(self.fill_auto_crop)
        crop_lay.addWidget(self.obj_crop, 1)
        crop_lay.addWidget(crop_btn)
        f.addRow('Crop (voxels)', self.crop_row)
        self.obj_pad = dspin(0.8, 0.0, 10.0, 0.1, 2)
        self.obj_pad.setToolTip('Medium around the object in the automatic crop (µm per side)')
        f.addRow('Crop margin (µm)', self.obj_pad)
        self.obj_radius = dspin(1.5, 0.1, 20, 0.1, 2)
        self.obj_dn = dspin(0.03, -1, 1, 0.005, 4)
        self.obj_voxel = dspin(0.1, 0.02, 0.5, 0.01, 3)
        f.addRow('Radius (µm)', self.obj_radius)
        f.addRow('Δn', self.obj_dn)
        f.addRow('Voxel (µm)', self.obj_voxel)
        self.n_medium = dspin(0, 0, 3, 0.001, 5, special='from object file')
        self.n_medium.setToolTip('0 = use the index stored in the object (.mat meta.n_medium, '
                                 '.h5 attribute n_medium). No index correction is applied.')
        f.addRow('n medium', self.n_medium)
        self.wavelength = dspin(0.6328, 0.2, 2.0, 0.001, 4)
        f.addRow('Wavelength (µm)', self.wavelength)
        self.obj_info = QLabel()
        self.obj_info.setObjectName('Info')
        self.obj_info.setWordWrap(True)
        f.addRow(self.obj_info)
        for w in (self.obj_radius, self.obj_dn, self.obj_voxel, self.n_medium, self.obj_pad):
            w.valueChanged.connect(lambda _: self.preview_timer.start())
        self.obj_file.changed.connect(self.preview_timer.start)
        self.obj_crop.editingFinished.connect(self.preview_timer.start)
        col.addWidget(g)

        # Simulation
        g = QGroupBox('Simulation')
        f = QFormLayout(g)
        self.resolution = ispin(12, 6, 40)
        self.polarization = QComboBox()
        self.polarization.addItems(['Ex', 'Ey'])
        self.pml = dspin(1.0, 0.2, 5, 0.1, 2)
        self.lat_margin = dspin(0.8, 0, 10, 0.1, 2)
        self.walkoff = dspin(1.0, 0, 10, 0.1, 2)
        self.z_margin = dspin(0.4, 0, 10, 0.1, 2)
        self.sensor_gap = dspin(0.4, 0.05, 10, 0.1, 2)
        self.mpi = ispin(4, 1, 64)
        self.max_mem = dspin(0, 0, 512, 1, 1, special='no limit')
        f.addRow('Resolution (px/µm)', self.resolution)
        f.addRow('Polarization', self.polarization)
        f.addRow('PML (µm)', self.pml)
        f.addRow('Lateral margin (µm)', self.lat_margin)
        f.addRow('Walk-off margin (µm)', self.walkoff)
        f.addRow('z margin (µm)', self.z_margin)
        f.addRow('Sensor gap (µm)', self.sensor_gap)
        f.addRow('MPI processes', self.mpi)
        f.addRow('Max. memory (GB)', self.max_mem)
        col.addWidget(g)

        # Illuminations
        g = QGroupBox('Illuminations')
        self.f_ill = f = QFormLayout(g)
        self.sel_mode = QComboBox()
        self.sel_mode.addItems(['Automatic (FPS, uniform)', 'From list file'])
        self.sel_mode.currentIndexChanged.connect(self.on_sel_mode)
        f.addRow('Selection', self.sel_mode)
        self.sel_n = ispin(32, 1, 2000)
        f.addRow('Number per view', self.sel_n)
        self.list_bottom = PathEdit('file', 'CSV (*.csv);;All files (*)')
        self.list_top = PathEdit('file', 'CSV (*.csv);;All files (*)')
        f.addRow('Bottom list', self.list_bottom)
        f.addRow('Top list', self.list_top)
        self.list_free = PathEdit('file', 'CSV (*.csv);;All files (*)',
                                  placeholder='inputs/illuminations/*.csv')
        f.addRow('List file', self.list_free)
        self.gen_btn = QPushButton('Generate list…')
        self.gen_btn.clicked.connect(self.generate_list)
        f.addRow(self.gen_btn)
        self.ill_info = QLabel()
        self.ill_info.setObjectName('Info')
        self.ill_info.setWordWrap(True)
        f.addRow(self.ill_info)
        for w in (self.list_bottom, self.list_top, self.list_free):
            w.changed.connect(self.update_illumination_info)
        col.addWidget(g)

        # Process
        g = QGroupBox('Process')
        v = QVBoxLayout(g)
        row = QHBoxLayout()
        self.do_trans = QCheckBox('Transmission')
        self.do_refl = QCheckBox('Reflection (4Pi)')
        self.do_trans.setChecked(True)
        for w in (self.do_trans, self.do_refl):
            w.toggled.connect(self.build_steps)
            row.addWidget(w)
        v.addLayout(row)
        self.steps = QListWidget()
        self.steps.setMinimumHeight(190)
        v.addWidget(self.steps)
        row = QHBoxLayout()
        row.addWidget(QLabel('Illuminations per view'))
        self.limit = ispin(0, 0, 10000, special='all')
        row.addWidget(self.limit)
        v.addLayout(row)
        self.warn = QLabel()
        self.warn.setObjectName('Warn')
        self.warn.setWordWrap(True)
        v.addWidget(self.warn)
        self.run_btn = QPushButton('▶  Run selected steps')
        self.run_btn.setObjectName('RunButton')
        self.run_btn.clicked.connect(self.run)
        v.addWidget(self.run_btn)
        self.stop_btn = QPushButton('■  Stop')
        self.stop_btn.setObjectName('StopButton')
        self.stop_btn.clicked.connect(self.runner.stop)
        self.stop_btn.setEnabled(False)
        v.addWidget(self.stop_btn)
        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.progress.setFormat('%p%')
        v.addWidget(self.progress)
        self.step_label = QLabel('Idle')
        self.step_label.setObjectName('Info')
        v.addWidget(self.step_label)
        b = QPushButton('📂  Open output folder')
        b.setObjectName('GreenButton')
        b.clicked.connect(self.open_output)
        v.addWidget(b)
        col.addWidget(g)
        col.addStretch()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setWidget(inner)
        scroll.setMinimumWidth(380)
        return scroll

    # ----------------------------------------------------------------- right panel
    def build_right(self):
        right = QSplitter(Qt.Vertical)
        self.tabs = QTabWidget()
        self.tabs.addTab(self.build_viewer(), 'Results viewer')
        self.tabs.addTab(self.build_figures(), 'Saved figures')
        right.addWidget(self.tabs)
        box = QWidget()
        lay = QVBoxLayout(box)
        lay.setContentsMargins(0, 4, 0, 0)
        t = QLabel('Console')
        t.setObjectName('PanelTitle')
        lay.addWidget(t)
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.console.setMaximumBlockCount(5000)
        lay.addWidget(self.console)
        right.addWidget(box)
        right.setSizes([800, 200])
        return right

    def build_viewer(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        bar = QHBoxLayout()
        t = QLabel('Results Viewer')
        t.setObjectName('AppTitle')
        bar.addWidget(t)
        bar.addStretch()
        self.v_view, self.v_name, self.v_source, self.v_cmap = (QComboBox() for _ in range(4))
        self.v_cmap.addItems(CMAPS)
        self.v_log = QCheckBox('Log scale')
        for label, combo in (('View', self.v_view), ('Illumination', self.v_name),
                             ('Hologram', self.v_source), ('Colormap', self.v_cmap)):
            bar.addWidget(QLabel(label + ':'))
            combo.setMinimumWidth(110)
            bar.addWidget(combo)
        bar.addWidget(self.v_log)
        b = QPushButton('⟳ Refresh')
        b.clicked.connect(self.refresh_viewer)
        bar.addWidget(b)
        lay.addLayout(bar)
        self.v_view.currentIndexChanged.connect(self.on_view)
        self.v_source.currentIndexChanged.connect(self.on_source_changed)
        self.v_name.currentIndexChanged.connect(self.on_name)
        self.v_cmap.currentIndexChanged.connect(self.on_name)
        self.v_log.toggled.connect(self.on_name)

        grid = QGridLayout()
        grid.setSpacing(8)
        self.p_obj = Panel('Material n(x, z)', tabs=('XZ', 'YZ', 'XY'))
        self.p_obj.tab_changed.connect(lambda _: self.update_object_panel(force=True))
        self.p_ill = Panel('Illuminations (kx/k, ky/k)')
        self.p_holo = Panel('Hologram', tabs=('Intensity', 'φ wrapped', 'φ unwrapped', 'Amplitude'))
        self.p_holo.tab_changed.connect(lambda _: self.draw_hologram())
        self.p_famp = Panel('Sensor field amplitude')
        self.p_fph = Panel('Sensor field phase')
        info = QFrame()
        info.setObjectName('Panel')
        il = QVBoxLayout(info)
        il.setContentsMargins(6, 6, 6, 6)
        t = QLabel('Details')
        t.setObjectName('PanelTitle')
        t.setAlignment(Qt.AlignCenter)
        il.addWidget(t)
        self.details = QTextBrowser()
        il.addWidget(self.details)
        for i, p in enumerate((self.p_obj, self.p_ill, self.p_holo, self.p_famp, self.p_fph, info)):
            grid.addWidget(p, i // 3, i % 3)
        lay.addLayout(grid, 1)
        return w

    def build_figures(self):
        w = QWidget()
        lay = QHBoxLayout(w)
        left = QVBoxLayout()
        self.fig_list = QListWidget()
        self.fig_list.currentItemChanged.connect(self.on_figure)
        self.fig_list.itemDoubleClicked.connect(
            lambda it: QDesktopServices.openUrl(QUrl.fromLocalFile(it.data(Qt.UserRole))))
        hint = QLabel('Figures saved by every step (PNG). Double-click: open with the system viewer.')
        hint.setObjectName('Muted')
        hint.setWordWrap(True)
        left.addWidget(hint)
        left.addWidget(self.fig_list)
        holder = QWidget()
        holder.setLayout(left)
        holder.setMaximumWidth(380)
        lay.addWidget(holder)
        self.p_fig = Panel('Figure', toolbar=True)
        lay.addWidget(self.p_fig, 1)
        return w

    # -------------------------------------------------------------- config <-> form
    def free(self):
        return self.workflow.currentIndex() == 1

    def load_config(self, path):
        try:
            cfg = json.loads(Path(path).read_text())
        except (OSError, ValueError) as e:
            QMessageBox.warning(self, 'Config', f'Could not read {path}:\n{e}')
            return
        self.cfg, self.cfg_path = cfg, Path(path)
        self._views_memory = dict(cfg.get('views', {}))
        self.cfg_to_form(cfg)
        self.cfg_label.setText(data.relative(path))
        self.log(f'Loaded config {data.relative(path)}')
        self.refresh_viewer()

    def cfg_to_form(self, cfg):
        o, s = cfg.get('object', {}), cfg.get('simulation', {})
        sel = cfg.get('selection', {})
        self.workflow.blockSignals(True)
        self.workflow.setCurrentIndex(1 if data.is_free(cfg) else 0)
        self.workflow.blockSignals(False)
        self.name.setText(cfg.get('name', ''))
        self.output.setText(cfg.get('output', ''))
        views = cfg.get('views', {})
        self.view_bottom.setText(views.get('bottom', ''))
        self.view_top.setText(views.get('top', ''))
        self.obj_source.setCurrentText(o.get('source', 'registered_gt'))
        self.obj_file.setText(o.get('file', ''))
        self.obj_crop.setText(json.dumps(o['crop_xyz']) if o.get('crop_xyz') else '')
        self.obj_pad.setValue(o.get('crop_pad_um', 0.8))
        self.obj_radius.setValue(o.get('radius_um', 1.5))
        self.obj_dn.setValue(o.get('dn', 0.03))
        self.obj_voxel.setValue(o.get('voxel_um', 0.1))
        self.n_medium.setValue(o.get('n_medium') or 0)
        self.wavelength.setValue(cfg.get('wavelength_um') or 0.6328)
        self.resolution.setValue(s.get('resolution', 12))
        self.polarization.setCurrentText(s.get('polarization', 'Ex'))
        self.pml.setValue(s.get('pml_um', 1.0))
        self.lat_margin.setValue(s.get('lateral_margin_um', 0.8))
        self.walkoff.setValue(s.get('walkoff_margin_um', 1.0))
        self.z_margin.setValue(s.get('z_margin_um', 0.4))
        self.sensor_gap.setValue(s.get('sensor_gap_um', 0.4))
        self.mpi.setValue(s.get('mpi_processes', 4))
        self.max_mem.setValue(s.get('max_memory_gb') or 0)
        lists = sel.get('lists', {})
        self.sel_mode.setCurrentIndex(1 if lists else 0)
        self.sel_n.setValue(sel.get('n', 32))
        self.list_bottom.setText(lists.get('bottom', ''))
        self.list_top.setText(lists.get('top', ''))
        self.list_free.setText(cfg.get('illuminations', ''))
        modes = cfg.get('modes', ['transmission'])
        if data.is_free(cfg):
            self.do_trans.setChecked('transmission' in modes)
            self.do_refl.setChecked('reflection' in modes)
        self.on_workflow()

    def form_to_cfg(self):
        cfg = copy.deepcopy(self.cfg)
        cfg['name'] = self.name.text().strip() or 'run'
        cfg['output'] = self.output.text() or f"outputs/{cfg['name']}"
        o = dict(cfg.get('object', {}))
        src = self.obj_source.currentText()
        for k in ('file', 'crop_xyz', 'crop_pad_um', 'n_immersion_gt', 'soft_mask_dn', 'radius_um',
                  'dn', 'voxel_um', 'edge_um', 'n_medium'):
            o.pop(k, None)
        o['source'] = src
        if src == 'sphere':
            o.update(radius_um=self.obj_radius.value(), dn=self.obj_dn.value(),
                     voxel_um=self.obj_voxel.value(), edge_um=self.cfg.get('object', {}).get('edge_um', 0.05))
        else:
            o['file'] = self.obj_file.text()
        if src == 'registered_gt':
            txt = self.obj_crop.text().strip()
            if txt:
                o['crop_xyz'] = json.loads(txt)
            else:
                o['crop_pad_um'] = self.obj_pad.value()
        if self.n_medium.value() > 0:
            o['n_medium'] = self.n_medium.value()
        cfg['object'] = o
        s = dict(cfg.get('simulation', {}))
        s.update(resolution=self.resolution.value(), polarization=self.polarization.currentText(),
                 pml_um=self.pml.value(), lateral_margin_um=self.lat_margin.value(),
                 walkoff_margin_um=self.walkoff.value(), z_margin_um=self.z_margin.value(),
                 sensor_gap_um=self.sensor_gap.value(), mpi_processes=self.mpi.value(),
                 max_memory_gb=self.max_mem.value() or None)
        cfg['simulation'] = s
        if self.free():
            cfg.pop('views', None)
            cfg['wavelength_um'] = self.wavelength.value()
            cfg['illuminations'] = self.list_free.text()
            cfg['modes'] = [m for m, c in (('transmission', self.do_trans), ('reflection', self.do_refl))
                            if c.isChecked()] or ['transmission']
        else:
            for k in ('wavelength_um', 'illuminations', 'modes'):
                cfg.pop(k, None)
            cfg['views'] = {v: w.text() for v, w in (('bottom', self.view_bottom), ('top', self.view_top))
                            if w.text()}
            sel = dict(cfg.get('selection', {}))
            sel['n'] = self.sel_n.value()
            lists = {v: w.text() for v, w in (('bottom', self.list_bottom), ('top', self.list_top))
                     if w.text()}
            sel['lists'] = lists if self.sel_mode.currentIndex() == 1 else {}
            cfg['selection'] = sel
        # Orden de claves estable y legible
        order = ['name', 'output', 'views', 'wavelength_um', 'illuminations', 'modes', 'object',
                 'simulation', 'selection', 'calibration', 'reflection', 'reference_views']
        return {k: cfg[k] for k in order if k in cfg} | {k: v for k, v in cfg.items() if k not in order}

    def write_config(self, path):
        cfg = self.form_to_cfg()
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(cfg, indent=2, ensure_ascii=False) + '\n')
        RunConfig.load(path)
        self.cfg, self.cfg_path = cfg, path
        self.cfg_label.setText(data.relative(path))
        self.statusBar().showMessage(f'Config saved: {data.relative(path)}', 5000)
        return path

    def open_config(self):
        p, _ = QFileDialog.getOpenFileName(self, 'Open config', str(ROOT / 'inputs' / 'config'),
                                           'JSON (*.json)')
        if p:
            self.load_config(p)

    def save_config(self):
        if not self.cfg_path:
            return self.save_config_as()
        try:
            return self.write_config(self.cfg_path)
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, 'Config', str(e))

    def save_config_as(self):
        start = ROOT / 'inputs' / 'config' / f"{self.name.text().strip() or 'run'}.json"
        p, _ = QFileDialog.getSaveFileName(self, 'Save config', str(start), 'JSON (*.json)')
        if p:
            try:
                return self.write_config(p)
            except (ValueError, OSError) as e:
                QMessageBox.warning(self, 'Config', str(e))

    # ------------------------------------------------------------- form behaviour
    def on_workflow(self):
        free = self.free()
        if free:
            self._views_memory = {'bottom': self.view_bottom.text(), 'top': self.view_top.text()}
        elif not self.view_bottom.text() and not self.view_top.text():
            self.view_bottom.setText(self._views_memory.get('bottom', ''))
            self.view_top.setText(self._views_memory.get('top', ''))
        self.g_data.setVisible(not free)
        set_row_visible(self.f_obj, self.wavelength, free)
        self.n_medium.setSpecialValueText('required' if free else 'from config_manip')
        self.do_refl.setText('Reflection' if free else 'Reflection (4Pi)')
        self.on_source()
        self.on_sel_mode()
        self.update_data_info()
        self.build_steps()

    def on_source(self):
        src = self.obj_source.currentText()
        for w in (self.obj_radius, self.obj_dn, self.obj_voxel):
            set_row_visible(self.f_obj, w, src == 'sphere')
        set_row_visible(self.f_obj, self.obj_file, src != 'sphere')
        for w in (self.crop_row, self.obj_pad):
            set_row_visible(self.f_obj, w, src == 'registered_gt')
        self.n_medium.setSpecialValueText('required' if src == 'sphere' and self.free()
                                          else 'from object file' if src != 'sphere'
                                          else 'from config_manip')
        self.preview_timer.start()

    def fill_auto_crop(self):
        cfg = self.viewer_cfg()
        cfg['object'].pop('crop_xyz', None)
        cfg['object']['crop_pad_um'] = self.obj_pad.value()
        try:
            self.statusBar().showMessage('Computing the object support…')
            QApplication.processEvents()
            crop = data.crop_for(cfg)
        except Exception as e:
            QMessageBox.warning(self, 'Auto crop', str(e))
            return
        self.obj_crop.setText(json.dumps(crop))
        self.statusBar().showMessage(f'Crop = {crop}', 5000)
        self.update_object_panel(force=True)

    def update_object_info(self, cfg, n=None, step=None):
        nm, origin = data.medium_info(cfg)
        lines = [f'n medium = {nm:.5f} ({origin})' if nm else '⚠ n medium not set']
        n0 = [i['n0'] for i in (data.acquisition_info(f) for f in cfg.get('views', {}).values()) if i]
        if nm and n0 and abs(n0[0] - nm) > 1e-4:
            lines.append(f'⚠ differs from N0 = {n0[0]:g} of config_manip (camera optics)')
        if n is not None:
            size = [s * m for s, m in zip(step, n.shape[::-1])]
            lines.append(f'box {n.shape[2]}×{n.shape[1]}×{n.shape[0]} vox = '
                         f'{size[0]:.1f}×{size[1]:.1f}×{size[2]:.1f} µm, centred at (0, 0, 0)')
            lines.append(f'n inside: {n.min():.4f} … {n.max():.4f}')
            edge = max(np.abs(n[[0, -1]] - nm).max(), np.abs(n[:, [0, -1]] - nm).max(),
                       np.abs(n[:, :, [0, -1]] - nm).max()) if nm else 0
            if edge > 1e-3:
                lines.append(f'⚠ object touches the box border (|n − n_medium| = {edge:.3g}): '
                             'enlarge the crop')
        self.obj_info.setText('\n'.join(lines))

    def on_sel_mode(self):
        free = self.free()
        from_list = self.sel_mode.currentIndex() == 1
        set_row_visible(self.f_ill, self.sel_mode, not free)
        set_row_visible(self.f_ill, self.sel_n, not free and not from_list)
        for w in (self.list_bottom, self.list_top):
            set_row_visible(self.f_ill, w, not free and from_list)
        set_row_visible(self.f_ill, self.list_free, free)
        self.gen_btn.setVisible(free)
        self.update_illumination_info()

    def update_data_info(self):
        lines = []
        for view, w in (('bottom', self.view_bottom), ('top', self.view_top)):
            info = data.acquisition_info(w.text()) if w.text() else None
            if info:
                lines.append(f"{view}: λ = {info['wavelength_um']:.4f} µm · N0 = {info['n0']} · "
                             f"NA = {info['NA']} · {info['dim']} px · {info['n_holograms']} holograms")
            elif w.text():
                lines.append(f'{view}: ⚠ no config_manip.txt in this folder')
        self.data_info.setText('\n'.join(lines) or 'Choose one folder per view.')

    def update_illumination_info(self):
        paths = [self.list_free.text()] if self.free() else (
            [self.list_bottom.text(), self.list_top.text()] if self.sel_mode.currentIndex() == 1 else [])
        lines = []
        for p in paths:
            if not p:
                continue
            try:
                rows = read_illumination_list(data.resolve(p))
                ks = [r for r in rows if 'kx' in r]
                th = max((r['theta_deg'] for r in ks), default=float('nan'))
                lines.append(f'{Path(p).name}: {len(rows)} illuminations'
                             + (f' · θmax = {th:.1f}° · ρmax = {np.sin(np.radians(th)):.3f}' if ks else
                                ' (by hologram file)'))
            except (OSError, ValueError) as e:
                lines.append(f'⚠ {Path(p).name}: {e}')
        if not self.free() and self.sel_mode.currentIndex() == 0:
            lines.append(f'The {self.sel_n.value()} most uniform directions (farthest-point '
                         'sampling) are taken from each view.')
        self.ill_info.setText('\n'.join(lines))

    def generate_list(self):
        dlg = GenerateDialog(self, f'inputs/illuminations/{self.name.text().strip() or "run"}_list.csv')
        if dlg.exec_() == QDialog.Accepted:
            self.list_free.setText(dlg.saved_path)
            self.update_illumination_info()
            self.refresh_viewer()

    def build_steps(self):
        self.steps.clear()
        t, r = self.do_trans.isChecked(), self.do_refl.isChecked()
        if self.free():
            items = [('prepare', 'Prepare object', True)]
            items += [('free:transmission', 'Transmission', True)] if t else []
            items += [('free:reflection', 'Reflection', True)] if r else []
            items += [('mie', 'Mie check of the solver (sphere, ~1 min)', False)]
            self.warn.setText('Small objects run in seconds. Memory and time grow with the '
                              'object volume and the resolution³.')
        else:
            items = [('prepare', 'Prepare object', True), ('select', 'Select illuminations', True),
                     ('estimate', 'Estimate memory per illumination', False),
                     ('calibrate', 'Calibrate object → camera', True)]
            items += [('run:bottom', 'Transmission bottom', True), ('run:top', 'Transmission top', True)] if t else []
            items += [('reflect:refl_top', 'Reflection refl_top (illum. bottom → camera top)', True),
                      ('reflect:refl_bottom', 'Reflection refl_bottom (illum. top → camera bottom)', True)] if r else []
            items += [('compare', 'Compare real vs simulated', True)] if t else []
            items += [('reference', 'Lateral reference views (XZ/YZ)', False)]
            self.warn.setText('Phantom at 12 px/µm: 2.5–5 min and 5–8 GB per illumination. '
                              'Close MATLAB and heavy apps. Steps resume: existing fields are skipped.')
        for key, label, checked in items:
            it = QListWidgetItem(label)
            it.setData(Qt.UserRole, key)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if checked else Qt.Unchecked)
            self.steps.addItem(it)

    # ------------------------------------------------------------------- running
    def step_args(self, key, cfg_path, cfg):
        c = ['--config', str(cfg_path)]
        lim = ['--limit', str(self.limit.value())] if self.limit.value() else []
        kind, _, arg = key.partition(':')
        if kind == 'run':
            return ['run'] + c + ['--view', arg] + lim
        if kind == 'reflect':
            return ['reflect'] + c + ['--view', arg] + lim
        if kind == 'free':
            return ['free'] + c + ['--mode', arg] + lim
        if kind == 'mie':
            return ['example-sphere', '--output', str(data.output_dir(cfg) / 'mie_check'),
                    '--resolution', str(self.resolution.value())]
        return [kind] + c

    def validate(self, cfg):
        problems = []
        if self.free():
            p = data.resolve(cfg.get('illuminations'))
            if not p or not p.exists():
                problems.append('Choose an illumination list file (or generate one).')
            if not data.n_medium_of(cfg):
                problems.append('Set n medium (the object does not store it and there is no real data).')
        else:
            if not cfg.get('views'):
                problems.append('Choose at least one view folder.')
            for v, folder in cfg.get('views', {}).items():
                if not data.acquisition_info(folder):
                    problems.append(f'View {v}: {folder} has no config_manip.txt.')
        if cfg['object']['source'] != 'sphere':
            p = data.resolve(cfg['object'].get('file'))
            if not p or not p.exists():
                problems.append(f"Object file not found: {cfg['object'].get('file')}")
        return problems

    def run(self):
        if self.runner.running():
            return
        keys = [self.steps.item(i).data(Qt.UserRole) for i in range(self.steps.count())
                if self.steps.item(i).checkState() == Qt.Checked]
        if not keys:
            QMessageBox.information(self, 'Run', 'Tick at least one step.')
            return
        try:
            cfg = self.form_to_cfg()
        except ValueError as e:
            QMessageBox.warning(self, 'Run', f'Invalid parameter: {e}')
            return
        problems = self.validate(cfg)
        if problems:
            QMessageBox.warning(self, 'Run', '\n'.join(problems))
            return
        path = self.cfg_path or ROOT / 'inputs' / 'config' / f"{cfg['name']}.json"
        try:
            path = self.write_config(path)
        except (ValueError, OSError) as e:
            QMessageBox.warning(self, 'Config', str(e))
            return
        labels = {self.steps.item(i).data(Qt.UserRole): self.steps.item(i).text()
                  for i in range(self.steps.count())}
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress.setValue(0)
        self.log(f'\n===== {cfg["name"]}: {len(keys)} step(s), config {data.relative(path)} =====')
        self.runner.run([(labels[k], self.step_args(k, path, cfg)) for k in keys])

    def on_progress(self, frac, label):
        self.progress.setValue(int(frac * 1000))
        self.step_label.setText(label)

    def on_finished(self, ok, msg):
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.step_label.setText(msg)
        self.statusBar().showMessage(msg)
        self.log(('✔ ' if ok else '✖ ') + msg)
        self.refresh_viewer()

    def log(self, text):
        self.console.appendPlainText(text)

    def open_output(self):
        out = data.output_dir(self.form_to_cfg())
        out.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(out)))

    # -------------------------------------------------------------------- viewer
    def viewer_cfg(self):
        try:
            return self.form_to_cfg()
        except ValueError:
            return self.cfg

    def refresh_viewer(self):
        cfg = self.viewer_cfg()
        self._holo_cache.clear()
        self.update_object_panel(force=True)
        self.update_figures(cfg)
        current = self.v_view.currentText()
        self.v_view.blockSignals(True)
        self.v_view.clear()
        self.v_view.addItems(data.available_views(cfg))
        if current:
            self.v_view.setCurrentText(current)
        self.v_view.blockSignals(False)
        self.on_view()

    def update_object_panel(self, force=False):
        cfg = self.viewer_cfg()
        key = json.dumps([cfg.get('object'), cfg.get('views'), str(data.output_dir(cfg)),
                          (data.output_dir(cfg) / 'object' / 'object.h5').exists()], default=str)
        if key != self._object_key or force:
            self._object_key = key
            try:
                self._object = data.load_object(cfg)
            except Exception as e:
                self._object = None
                self.p_obj.message(f'Object not available:\n{e}')
                self.update_object_info(cfg)
                return
        if not self._object:
            return
        n, step, label = self._object
        self.update_object_info(cfg, n, step)
        plane = ('XZ', 'YZ', 'XY')[self.p_obj.tab()]
        img, ext, xl, yl = data.object_slice(n, step, plane)
        self.p_obj.set_title(f'Material n({plane[0].lower()}, {plane[1].lower()})')
        ax = self.p_obj.image(img, ext, cmap='gray', xlabel=xl, ylabel=yl, clabel='n')
        ax.set_title(f'{label} · {n.shape[2]}×{n.shape[1]}×{n.shape[0]} vox', fontsize=7)
        self.p_obj.canvas.draw_idle()

    def illumination_key(self, cfg, view):
        if data.is_free(cfg):
            return 'list'
        refl = cfg.get('reflection', {}).get('views', {'refl_top': {'illumination': 'bottom'},
                                                        'refl_bottom': {'illumination': 'top'}})
        return refl[view]['illumination'] if view in refl else view

    def update_illumination_panel(self, cfg, view, name):
        try:
            sets = data.illumination_sets(cfg)
        except (OSError, ValueError) as e:
            self.p_ill.message(f'Illumination list error:\n{e}')
            return
        if not sets:
            self.p_ill.message('No illuminations yet.\nLoad a list file, or run "Select".')
            return
        self.p_ill.scatter_sets(sets, highlight=(self.illumination_key(cfg, view), name) if view else None)

    def on_view(self):
        cfg = self.viewer_cfg()
        view = self.v_view.currentText()
        srcs = data.hologram_sources(cfg, view) if view else {}
        self._sources = srcs
        current = self.v_source.currentText()
        self.v_source.blockSignals(True)
        self.v_source.clear()
        self.v_source.addItems(list(srcs))
        if current in srcs:
            self.v_source.setCurrentText(current)
        elif 'sim' in srcs:
            self.v_source.setCurrentText('sim')
        self.v_source.blockSignals(False)
        self.on_source_changed()

    def on_source_changed(self):
        cfg = self.viewer_cfg()
        view = self.v_view.currentText()
        folder = self._sources.get(self.v_source.currentText()) if view else None
        names = data.illumination_names(cfg, view, folder) if view else []
        current = self.v_name.currentText()
        self.v_name.blockSignals(True)
        self.v_name.clear()
        self.v_name.addItems(names)
        if current in names:
            self.v_name.setCurrentText(current)
        self.v_name.blockSignals(False)
        self.on_name()

    def on_name(self, *_):
        cfg = self.viewer_cfg()
        view, name = self.v_view.currentText(), self.v_name.currentText()
        self.update_illumination_panel(cfg, view, name)
        self.draw_hologram()
        self.draw_field(cfg, view, name)
        self.update_details(cfg, view, name)

    def draw_hologram(self):
        name = self.v_name.currentText()
        folder = self._sources.get(self.v_source.currentText()) if self.v_view.currentText() else None
        if not folder or not name or not (folder / f'{name}.pgm').exists():
            self.p_holo.message('Free mode: no camera holograms.\nSee the sensor field panels.'
                                if self.free() else 'No hologram for this view / illumination.')
            return
        key = (str(folder), name)
        if key not in self._holo_cache:
            try:
                cal = self.cfg.get('calibration', {})
                self._holo_cache[key] = data.hologram_maps(
                    folder, name, self.cfg.get('selection', {}).get('tukey', 0.1),
                    tuple(cal.get('background_annulus_um', (17, 33))), cal.get('background_poly_order', 4))
            except Exception as e:
                self.p_holo.message(f'Cannot read hologram:\n{e}')
                return
        m = self._holo_cache[key]
        tab = self.p_holo.tab()
        cmap = self.v_cmap.currentText()
        title = f'{self.v_source.currentText()} · {name}.pgm'
        if tab == 0:
            ax = self.p_holo.image(*m['intensity'], cmap=cmap, clabel='gray level', log=self.v_log.isChecked())
        elif tab == 1:
            ax = self.p_holo.image(*m['phase'], cmap='twilight', clabel='rad', vmin=-np.pi, vmax=np.pi)
        elif tab == 2:
            ax = self.p_holo.image(*m['unwrapped'], cmap=cmap, clabel='rad')
        else:
            ax = self.p_holo.image(*m['amplitude'], cmap=cmap, clabel='relative', log=self.v_log.isChecked())
        pk = m['peak']
        ax.set_title(f"{title} · θ = {pk['theta_deg']:.1f}°, φ = {pk['phi_deg']:.0f}°", fontsize=7)
        self.p_holo.canvas.draw_idle()

    def draw_field(self, cfg, view, name):
        npz = data.output_dir(cfg) / view / 'fields' / f'{name}.npz' if view and name else None
        if not npz or not npz.exists():
            msg = 'No sensor field yet.\nRun the simulation for this view.'
            self.p_famp.message(msg)
            self.p_fph.message(msg)
            return
        try:
            f = data.field_maps(npz)
        except Exception as e:
            self.p_famp.message(f'Cannot read field:\n{e}')
            self.p_fph.message('')
            return
        cmap = self.v_cmap.currentText()
        ax = self.p_famp.image(f['amplitude'], f['extent'], cmap=cmap, clabel=f['label'],
                               log=self.v_log.isChecked())
        ax.set_title(f"{view} {name} · {f['label']} ({f['component']}) at z = {f['z']:+.2f} µm", fontsize=7)
        ax = self.p_fph.image(f['phase'], f['extent'], cmap='twilight', clabel='rad', vmin=-np.pi, vmax=np.pi)
        ax.set_title(f"phase of {f['label'].strip('|')}", fontsize=7)
        self.p_famp.canvas.draw_idle()
        self.p_fph.canvas.draw_idle()

    def update_details(self, cfg, view, name):
        rows = []
        if view and name:
            rec = data.report_rows(cfg, view).get(name)
            if rec:
                keys = ('original', 'filename', 'theta_deg', 'phi_deg', 'theta_inc_deg', 'kx', 'ky',
                        'kx_norm_raw', 'ky_norm_raw', 'kx_norm_obj', 'ky_norm_obj', 'delta_p_bins',
                        'peak_to_specular_bins', 'rms_S', 'max_S', 'back_rms', 'saturated_fraction',
                        'n_cells_M', 'meep_time', 'wall_s')
                for k in keys:
                    if k in rec and rec[k] not in ('', None):
                        v = rec[k]
                        try:
                            v = f'{float(v):.4g}'
                        except ValueError:
                            pass
                        rows.append((k, v))
        html = [f'<h3 style="color:#ff6a2b">{view or "—"} · {name or "—"}</h3>']
        if rows:
            html.append('<table cellspacing="2">' + ''.join(
                f'<tr><td style="color:#8d93a0">{k}</td><td>&nbsp;{v}</td></tr>' for k, v in rows) + '</table>')
        else:
            html.append('<p style="color:#8d93a0">No report for this illumination yet.</p>')
        cal = data.calibration_summary(cfg)
        if cal:
            html.append('<h4 style="color:#ff6a2b">Calibration</h4>')
            for v, b in cal.items():
                html.append(f"<p>{v}: entry {b['orientation']}, M = {b['M']}, conj = {b['conj']}, "
                            f"corr = {b['corr']:.3f}, focus = {b['z_focus_um']:+.2f} µm</p>")
        html.append(f'<p style="color:#8d93a0">Output: {data.relative(data.output_dir(cfg))}</p>')
        self.details.setHtml(''.join(html))

    def update_figures(self, cfg):
        current = self.fig_list.currentItem().data(Qt.UserRole) if self.fig_list.currentItem() else None
        self.fig_list.blockSignals(True)
        self.fig_list.clear()
        out = data.output_dir(cfg)
        for p in data.saved_figures(cfg):
            it = QListWidgetItem(str(p.relative_to(out)))
            it.setData(Qt.UserRole, str(p))
            self.fig_list.addItem(it)
            if str(p) == current:
                self.fig_list.setCurrentItem(it)
        self.fig_list.blockSignals(False)
        if self.fig_list.count() and not self.fig_list.currentItem():
            self.fig_list.setCurrentRow(0)
        elif not self.fig_list.count():
            self.p_fig.message('No saved figures in this output folder yet.')

    def on_figure(self, item, _prev=None):
        if item:
            self.p_fig.set_title(item.text())
            self.p_fig.picture(item.data(Qt.UserRole))

    def closeEvent(self, event):
        if self.runner.running():
            if QMessageBox.question(self, 'Quit', 'A simulation is running. Stop it and quit?') != QMessageBox.Yes:
                event.ignore()
                return
            self.runner.stop()
        event.accept()
