# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Dark theme (Qt style sheet + matplotlib) with an orange accent.

BG = '#15171c'
PANEL = '#1d2027'
FIELD = '#262a33'
BORDER = '#343a46'
TEXT = '#d8dbe1'
MUTED = '#8d93a0'
ACCENT = '#ff6a2b'
ACCENT_DARK = '#c94e1a'
GREEN = '#3ecf6e'

QSS = f"""
QWidget {{ background: {BG}; color: {TEXT}; font-size: 13px; }}
QScrollArea, QScrollArea > QWidget > QWidget {{ background: {BG}; border: none; }}
QGroupBox {{ background: {PANEL}; border: 1px solid {BORDER}; border-radius: 8px;
             margin-top: 14px; padding: 12px 8px 8px 8px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 12px; padding: 0 4px;
                    color: {ACCENT}; font-weight: bold; font-size: 14px; }}
QGroupBox QWidget {{ background: {PANEL}; }}
QLabel#AppTitle {{ font-size: 18px; font-weight: bold; color: {TEXT}; }}
QLabel#Version {{ color: {ACCENT}; font-weight: bold; }}
QLabel#Muted, QLabel#Info {{ color: {MUTED}; font-size: 12px; }}
QLabel#PanelTitle {{ color: {ACCENT}; font-weight: bold; font-size: 14px; padding: 4px;
                     background: {FIELD}; border-radius: 4px; }}
QLabel#Warn {{ color: #ffb347; font-size: 12px; }}
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{
    background: {FIELD}; border: 1px solid {BORDER}; border-radius: 5px; padding: 4px 6px;
    selection-background-color: {ACCENT}; min-height: 20px; }}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{ border: 1px solid {ACCENT}; }}
QLineEdit:disabled, QSpinBox:disabled, QDoubleSpinBox:disabled, QComboBox:disabled {{ color: {MUTED}; }}
QComboBox QAbstractItemView {{ background: {FIELD}; border: 1px solid {BORDER};
                               selection-background-color: {ACCENT}; }}
QPushButton {{ background: {FIELD}; border: 1px solid {BORDER}; border-radius: 6px;
               padding: 6px 10px; }}
QPushButton:hover {{ border: 1px solid {ACCENT}; }}
QPushButton:pressed {{ background: {BORDER}; }}
QPushButton:disabled {{ color: {MUTED}; border-color: {FIELD}; }}
QPushButton#RunButton {{ border: 2px solid {ACCENT}; font-size: 15px; font-weight: bold;
                         padding: 10px; }}
QPushButton#RunButton:hover {{ background: {ACCENT_DARK}; }}
QPushButton#StopButton {{ border: 2px solid #d9534f; font-weight: bold; padding: 10px; }}
QPushButton#GreenButton {{ border: 2px solid {GREEN}; color: {GREEN}; font-weight: bold; }}
QProgressBar {{ background: {FIELD}; border: 1px solid {BORDER}; border-radius: 5px;
                text-align: center; height: 18px; }}
QProgressBar::chunk {{ background: {ACCENT}; border-radius: 4px; }}
QCheckBox {{ spacing: 6px; }}
QCheckBox::indicator {{ width: 15px; height: 15px; border: 1px solid {BORDER}; border-radius: 3px;
                        background: {FIELD}; }}
QCheckBox::indicator:checked {{ background: {ACCENT}; border-color: {ACCENT}; }}
QListWidget {{ background: {FIELD}; border: 1px solid {BORDER}; border-radius: 5px; }}
QListWidget::item {{ padding: 3px; }}
QListWidget::item:selected {{ background: {ACCENT}; color: white; }}
QTabWidget::pane {{ border: 1px solid {BORDER}; border-radius: 6px; background: {PANEL}; }}
QTabBar::tab {{ background: {FIELD}; color: {MUTED}; padding: 6px 14px; border: 1px solid {BORDER};
                border-bottom: none; border-top-left-radius: 5px; border-top-right-radius: 5px;
                margin-right: 2px; }}
QTabBar::tab:selected {{ background: {ACCENT}; color: white; }}
QFrame#Panel {{ background: {PANEL}; border: 1px solid {BORDER}; border-radius: 8px; }}
QFrame#Panel QWidget {{ background: {PANEL}; }}
QPlainTextEdit, QTextBrowser {{ background: #101216; border: 1px solid {BORDER}; border-radius: 6px;
                                font-family: Menlo, Monaco, monospace; font-size: 12px; }}
QSplitter::handle {{ background: {BORDER}; }}
QScrollBar:vertical {{ background: {BG}; width: 10px; }}
QScrollBar::handle:vertical {{ background: {BORDER}; border-radius: 5px; min-height: 24px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; }}
QStatusBar {{ color: {MUTED}; }}
QToolTip {{ background: {FIELD}; color: {TEXT}; border: 1px solid {ACCENT}; }}
"""


def apply_matplotlib_style():
    import matplotlib
    matplotlib.rcParams.update({
        'figure.facecolor': PANEL, 'axes.facecolor': PANEL, 'savefig.facecolor': PANEL,
        'axes.edgecolor': MUTED, 'axes.labelcolor': TEXT, 'text.color': TEXT,
        'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.titlecolor': TEXT,
        'font.size': 9, 'axes.titlesize': 9, 'axes.labelsize': 9,
        'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
        'legend.facecolor': FIELD, 'legend.edgecolor': BORDER, 'grid.color': BORDER,
        'image.cmap': 'hot',
    })
