# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Matplotlib panels embedded in Qt (title bar, optional tabs, colorbar).

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import QFrame, QLabel, QTabBar, QVBoxLayout

from . import style


class Panel(QFrame):
    tab_changed = pyqtSignal(int)

    def __init__(self, title, tabs=(), toolbar=False, parent=None):
        super().__init__(parent)
        self.setObjectName('Panel')
        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(4)
        self.title = QLabel(title)
        self.title.setObjectName('PanelTitle')
        self.title.setAlignment(Qt.AlignCenter)
        lay.addWidget(self.title)
        self.tabs = None
        if tabs:
            self.tabs = QTabBar()
            self.tabs.setExpanding(False)
            self.tabs.setDrawBase(False)
            for t in tabs:
                self.tabs.addTab(t)
            self.tabs.currentChanged.connect(self.tab_changed)
            lay.addWidget(self.tabs)
        self.fig = Figure(figsize=(4, 3.4), layout='constrained')
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setMinimumSize(220, 200)
        if toolbar:
            lay.addWidget(NavigationToolbar2QT(self.canvas, self))
        lay.addWidget(self.canvas, 1)
        self.message('—')

    def tab(self):
        return self.tabs.currentIndex() if self.tabs else 0

    def set_title(self, text):
        self.title.setText(text)

    def _axes(self):
        self.fig.clear()
        return self.fig.add_subplot(111)

    def message(self, text):
        ax = self._axes()
        ax.axis('off')
        ax.text(0.5, 0.5, text, ha='center', va='center', color=style.MUTED, wrap=True,
                transform=ax.transAxes)
        self.canvas.draw_idle()

    def image(self, img, extent=None, cmap='hot', xlabel='x (µm)', ylabel='y (µm)', clabel='',
              vmin=None, vmax=None, log=False):
        ax = self._axes()
        data = np.asarray(img, float)
        if log:
            data = np.log10(np.abs(data) + 1e-12)
        im = ax.imshow(data, origin='lower', extent=extent, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation='nearest', aspect='equal')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cb = self.fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.ax.tick_params(labelsize=7)
        if clabel:
            cb.set_label(clabel, fontsize=8)
        self.canvas.draw_idle()
        return ax

    def scatter_sets(self, sets, highlight=None, rho_max=None):
        ax = self._axes()
        colors = [style.ACCENT, '#4fc3f7', '#b388ff', style.GREEN]
        ci = 0
        for label, kx, ky, names, kind in sets:
            if kind == 'all':
                ax.plot(kx, ky, '.', color='#5b616d', ms=3, label=label)
                continue
            c = colors[ci % len(colors)]
            ci += 1
            ax.plot(kx, ky, 'o', color=c, ms=4.5, mfc='none' if ci > 1 else c, label=label)
            if len(names) <= 40:
                for x, y, n in zip(kx, ky, names):
                    ax.annotate(n.replace('i0', '').lstrip('0') or '0', (x, y), fontsize=6,
                                color=style.MUTED, xytext=(2, 2), textcoords='offset points')
            if highlight and highlight[0] == label.split(':')[0] and highlight[1] in names:
                j = names.index(highlight[1])
                ax.plot(kx[j], ky[j], 'o', ms=12, mfc='none', mec='white', mew=1.5)
        t = np.linspace(0, 2 * np.pi, 200)
        for r in (rho_max,) if rho_max else ():
            ax.plot(r * np.cos(t), r * np.sin(t), ':', color=style.MUTED, lw=0.8)
        ax.plot(np.cos(t), np.sin(t), '-', color=style.BORDER, lw=0.8)
        ax.set_aspect('equal')
        ax.set_xlim(-1.02, 1.02)
        ax.set_ylim(-1.02, 1.02)
        ax.grid(alpha=0.25)
        ax.set_xlabel('kx / k')
        ax.set_ylabel('ky / k')
        if sets:
            ax.legend(loc='lower left', fontsize=7, markerscale=0.8)
        self.canvas.draw_idle()

    def picture(self, path):
        import matplotlib.image as mpimg
        ax = self._axes()
        ax.imshow(mpimg.imread(str(path)))
        ax.axis('off')
        self.canvas.draw_idle()
