# FDTD_meep v1.1.0 | 2026-08-27 | Carlos Chacón
# Command line. v1.2: adds 'free' (any object + illumination list, no real data).

import argparse
import csv
import sys

from .config import RunConfig, ROOT


def main(argv=None):
    parser = argparse.ArgumentParser(prog='python -m fdtd_meep',
                                     description='Simulación FDTD (Meep, campo dispersado) de '
                                                 'adquisiciones MTD y hologramas sintéticos.')
    sub = parser.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('example-sphere', help='Esfera frente a Mie exacto (sin datos).')
    p.add_argument('--output', default=str(ROOT / 'outputs' / 'example_sphere'))
    p.add_argument('--resolution', type=int, default=12)
    p.add_argument('--reflection', action='store_true',
                   help='Retrodispersión (sensor bajo la esfera) frente a Mie.')

    for name, help_ in (('prepare', 'Prepara el objeto (outputs/<run>/object/object.h5).'),
                        ('select', 'Picos especulares, selección FPS y real_subset/.'),
                        ('estimate', 'Celdas y memoria estimadas por iluminación.'),
                        ('calibrate', 'Calibración objeto -> cámara con los hologramas reales.'),
                        ('compare', 'Campos demodulados reales frente a simulados.')):
        q = sub.add_parser(name, help=help_)
        q.add_argument('--config', required=True)

    q = sub.add_parser('run', help='Simula las iluminaciones de una vista y escribe sim_subset/.')
    q.add_argument('--config', required=True)
    q.add_argument('--view', required=True)
    q.add_argument('--limit', type=int)

    q = sub.add_parser('reflect', help='Reflexión: retrodispersión registrada por la cámara '
                                       'opuesta, con las iluminaciones de otra vista.')
    q.add_argument('--config', required=True)
    q.add_argument('--view', required=True, help='Vista de reflexión (reflection.views).')
    q.add_argument('--limit', type=int)

    q = sub.add_parser('free', help='Modo libre: objeto + lista de iluminaciones, sin hologramas '
                                    'reales (outputs/<run>/<modo>/).')
    q.add_argument('--config', required=True)
    q.add_argument('--mode', choices=('transmission', 'reflection'),
                   help='Por defecto, los "modes" de la configuración.')
    q.add_argument('--limit', type=int)

    q = sub.add_parser('reference', help='Figuras laterales de referencia (campo en el dominio FDTD, '
                                         'iluminación central). Solo se ejecuta si se pide.')
    q.add_argument('--config', required=True)
    q.add_argument('--view', action='append', help='Vista (repetible); por defecto todas.')

    q = sub.add_parser('sim', help='Una iluminación (uso interno, compatible con mpirun).')
    q.add_argument('--config', required=True)
    q.add_argument('--orientation', choices=('zmin', 'zmax'), required=True)
    q.add_argument('--k', required=True)
    q.add_argument('--out', required=True)
    q.add_argument('--reflection', action='store_true')
    q.add_argument('--slices', default='', help='Cortes laterales a guardar: XZ, YZ o XZ,YZ.')
    q.add_argument('--z-pad', type=float, default=0.0, help='Margen extra en z (µm) para los cortes.')
    q.add_argument('--slice-center', default=None, help='x,y,z del corte en el marco simulado (µm).')

    args = parser.parse_args(argv)
    if args.cmd == 'example-sphere':
        from .example_sphere import run
        output = args.output
        if args.reflection and output == str(ROOT / 'outputs' / 'example_sphere'):
            output += '_reflection'
        run(output, args.resolution, reflection=args.reflection)
        return

    from .pipeline import Pipeline
    pipe = Pipeline(RunConfig.load(args.config), args.config)
    if args.cmd == 'prepare':
        print(pipe.prepare())
    elif args.cmd == 'select':
        for view, chosen in pipe.select().items():
            print(f'{view}: {len(chosen)} iluminaciones -> {pipe.view_dir(view)}')
    elif args.cmd == 'estimate':
        rows = pipe.estimate()
        w = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    elif args.cmd == 'calibrate':
        report = pipe.calibrate()
        for view in pipe.acq:
            b = report[view]['best']
            print(f"{view}: entrada {b['orientation']}, M={b['M']}, conj={b['conj']}, "
                  f"t={b['t_um']}, foco={b['z_focus_um']:+.2f} µm")
    elif args.cmd == 'run':
        pipe.run_view(args.view, args.limit)
        if pipe.cfg.reference_views.enabled:
            from .reference_views import run_reference
            run_reference(pipe, [args.view])
    elif args.cmd == 'reflect':
        from .reflection import run_reflection
        run_reflection(pipe, args.view, args.limit)
        if pipe.cfg.reference_views.enabled:
            from .reference_views import run_reference
            run_reference(pipe, [args.view])
    elif args.cmd == 'free':
        from .free import run_free
        for mode in [args.mode] if args.mode else pipe.cfg.modes:
            records = run_free(pipe, mode, args.limit)
            print(f'{mode}: {len(records)} iluminaciones -> {pipe.view_dir(mode)}')
    elif args.cmd == 'reference':
        from .reference_views import run_reference
        for path in run_reference(pipe, args.view):
            print(path)
    elif args.cmd == 'compare':
        for r in pipe.compare():
            print(r)
    elif args.cmd == 'sim':
        kx, ky = map(float, args.k.split(','))
        pipe.simulate_one(args.orientation, kx, ky, args.out, args.reflection,
                          tuple(p for p in args.slices.split(',') if p), args.z_pad,
                          None if args.slice_center is None else [float(v) for v in args.slice_center.split(',')])


if __name__ == '__main__':
    main()
