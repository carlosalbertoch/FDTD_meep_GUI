# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Launch the GUI:  python -m fdtd_meep_gui [--config inputs/config/x.json]

import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(prog='python -m fdtd_meep_gui',
                                     description='FDTD_meep graphical interface (PyQt5).')
    parser.add_argument('--config', help='JSON config to open (default: the sphere example).')
    parser.add_argument('--screenshot', help='Save a PNG of the window and quit (docs/tests).')
    parser.add_argument('--tab', type=int, default=0, help='Right-hand tab for --screenshot.')
    args = parser.parse_args(argv)

    import matplotlib
    matplotlib.use('Qt5Agg')
    from PyQt5.QtCore import QLocale, QTimer
    from PyQt5.QtWidgets import QApplication

    from .app import MainWindow
    from .style import QSS, apply_matplotlib_style

    # Punto decimal en todos los campos, como en los JSON de configuración
    QLocale.setDefault(QLocale(QLocale.English, QLocale.UnitedStates))
    app = QApplication(sys.argv[:1])
    app.setApplicationName('FDTD_meep GUI')
    app.setStyle('Fusion')
    app.setStyleSheet(QSS)
    apply_matplotlib_style()
    win = MainWindow(args.config)
    win.show()
    if args.screenshot:
        win.tabs.setCurrentIndex(args.tab)

        def shoot():
            win.grab().save(args.screenshot)
            app.quit()
        QTimer.singleShot(2500, shoot)
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main())
