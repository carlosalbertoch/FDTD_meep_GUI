# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Run configuration (JSON). v1.2: sphere objects, illumination lists, free mode.

import json
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[1]


def _resolve(path):
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def _build(cls, data):
    names = {f.name for f in fields(cls)}
    unknown = set(data) - names
    if unknown:
        raise ValueError(f'{cls.__name__}: parámetros desconocidos {sorted(unknown)}')
    return cls(**data)


@dataclass
class ObjectParams:
    source: str = 'registered_gt'
    file: str = 'inputs/phantom/phantom_ground_truth_registered.mat'
    crop_xyz: Optional[List[Tuple[int, int]]] = None   # None = soporte + crop_pad_um
    crop_pad_um: float = 0.8
    n_medium: Optional[float] = None                   # None = el que indica el objeto
    n_immersion_gt: Optional[float] = None             # solo si se pide: corrección v1.1
    soft_mask_dn: float = 0.01
    center_um: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    # source = 'sphere'
    radius_um: float = 1.5
    dn: float = 0.03
    voxel_um: float = 0.1
    edge_um: float = 0.05


@dataclass
class SimulationParams:
    resolution: int = 12
    pml_um: float = 1.0
    pml_far_um: Optional[float] = None
    lateral_margin_um: float = 0.8
    walkoff_margin_um: float = 1.0
    z_margin_um: float = 0.4
    sensor_gap_um: float = 0.4
    pulse_fwidth: float = 0.4
    field_decay_by: float = 1e-8
    max_runtime_periods: float = 500.0
    polarization: str = 'Ex'
    mpi_processes: int = 4
    max_memory_gb: Optional[float] = None
    bytes_per_cell: float = 139.0


@dataclass
class SelectionParams:
    n: int = 32
    tukey: float = 0.10
    lists: Dict[str, str] = field(default_factory=dict)


@dataclass
class CalibrationParams:
    z_range_um: Tuple[float, float] = (-14.0, 14.0)
    z_step_um: float = 1.0
    z_refine_um: float = 0.25
    background_annulus_um: Tuple[float, float] = (17.0, 33.0)
    background_poly_order: int = 4
    match_radius_um: float = 18.0


@dataclass
class ReflectionParams:
    views: Dict[str, Dict[str, str]] = field(default_factory=lambda: {
        'refl_top': {'illumination': 'bottom', 'camera': 'top'},
        'refl_bottom': {'illumination': 'top', 'camera': 'bottom'}})
    reference_ratio: float = 2.0
    parasitic_specular: List[float] = field(default_factory=list)


@dataclass
class ReferenceViewParams:
    enabled: bool = False
    planes: List[str] = field(default_factory=lambda: ['YZ'])
    position_um: Optional[float] = None
    field_modes: List[str] = field(default_factory=lambda: ['total', 'scattered'])
    intensity_mode: str = 'compressed_arcsinh'
    arcsinh_alpha: Optional[float] = None
    colormap: str = 'turbo'
    overlay_contour: bool = True
    save_png: bool = True
    save_svg: bool = True
    z_pad_um: float = 3.0


@dataclass
class RunConfig:
    name: str
    output: Path
    views: Dict[str, Path]
    object: ObjectParams = field(default_factory=ObjectParams)
    simulation: SimulationParams = field(default_factory=SimulationParams)
    selection: SelectionParams = field(default_factory=SelectionParams)
    calibration: CalibrationParams = field(default_factory=CalibrationParams)
    reflection: ReflectionParams = field(default_factory=ReflectionParams)
    reference_views: ReferenceViewParams = field(default_factory=ReferenceViewParams)
    # Modo libre (sin hologramas reales): λ, lista de iluminaciones y modos a simular
    wavelength_um: Optional[float] = None
    illuminations: Optional[str] = None
    modes: List[str] = field(default_factory=lambda: ['transmission'])

    @classmethod
    def load(cls, path):
        data = json.loads(Path(path).read_text())
        name = data.get('name', Path(path).stem)
        return cls(
            name=name,
            output=_resolve(data.get('output', f'outputs/{name}')),
            views={k: _resolve(v) for k, v in data.get('views', {}).items()},
            wavelength_um=data.get('wavelength_um'),
            illuminations=data.get('illuminations'),
            modes=data.get('modes', ['transmission']),
            object=_build(ObjectParams, data.get('object', {})),
            simulation=_build(SimulationParams, data.get('simulation', {})),
            selection=_build(SelectionParams, data.get('selection', {})),
            calibration=_build(CalibrationParams, data.get('calibration', {})),
            reflection=_build(ReflectionParams, data.get('reflection', {})),
            reference_views=_build(ReferenceViewParams, data.get('reference_views', {})),
        )

    def object_file(self):
        return _resolve(self.object.file)

    def illuminations_file(self):
        if not self.illuminations:
            raise ValueError('Falta "illuminations": la lista de iluminaciones (CSV).')
        return _resolve(self.illuminations)

    def selection_list(self, view):
        path = self.selection.lists.get(view)
        return _resolve(path) if path else None
