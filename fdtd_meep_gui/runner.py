# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Runs a queue of `python -m fdtd_meep ...` steps in a separate process group (MPI-safe stop).

import os
import re
import signal
import sys

from PyQt5.QtCore import QObject, QProcess, QProcessEnvironment, QTimer, pyqtSignal

from fdtd_meep.config import ROOT

# New session so that Stop also ends mpirun and its ranks.
LAUNCH = ("import os, sys; os.setsid(); "
          "os.execv(sys.executable, [sys.executable, '-u', '-m', 'fdtd_meep'] + sys.argv[1:])")
PROGRESS = re.compile(r'^\[(\S+) (\d+)/(\d+)\]')


class Runner(QObject):
    line = pyqtSignal(str)
    step_started = pyqtSignal(int, int, str)
    progress = pyqtSignal(float, str)
    step_done = pyqtSignal(str)
    finished = pyqtSignal(bool, str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.proc = None
        self.steps = []
        self.index = 0
        self.stopped = False
        self._buffer = ''

    def running(self):
        return self.proc is not None

    def run(self, steps):
        """steps: [(etiqueta, [argumentos de fdtd_meep])]."""
        self.steps, self.index, self.stopped = list(steps), 0, False
        self._next()

    def _next(self):
        if self.index >= len(self.steps):
            self.proc = None
            self.progress.emit(1.0, 'Done')
            self.finished.emit(True, f'{len(self.steps)} step(s) completed.')
            return
        label, args = self.steps[self.index]
        self.step_started.emit(self.index, len(self.steps), label)
        self.progress.emit(self.index / len(self.steps), label)
        self.line.emit(f'\n▶ {label}\n$ python -m fdtd_meep {" ".join(args)}')
        p = QProcess(self)
        env = QProcessEnvironment.systemEnvironment()
        env.insert('PYTHONUNBUFFERED', '1')
        env.insert('MPLBACKEND', 'Agg')
        p.setProcessEnvironment(env)
        p.setWorkingDirectory(str(ROOT))
        p.setProcessChannelMode(QProcess.MergedChannels)
        p.readyReadStandardOutput.connect(self._read)
        p.finished.connect(self._finished)
        self.proc = p
        p.start(sys.executable, ['-c', LAUNCH] + list(args))

    def _read(self):
        self._buffer += bytes(self.proc.readAllStandardOutput()).decode('utf-8', 'replace')
        *lines, self._buffer = self._buffer.split('\n')
        for ln in lines:
            self.line.emit(ln)
            m = PROGRESS.match(ln)
            if m:
                i, n = int(m.group(2)), int(m.group(3))
                frac = (self.index + (i - 1) / n) / len(self.steps)
                self.progress.emit(frac, f'{self.steps[self.index][0]}: {i}/{n}')

    def _finished(self, code, status):
        if self._buffer:
            self.line.emit(self._buffer)
            self._buffer = ''
        label = self.steps[self.index][0]
        self.proc = None
        if self.stopped:
            self.finished.emit(False, f'Stopped during: {label}')
        elif code != 0 or status != QProcess.NormalExit:
            self.finished.emit(False, f'{label} failed (exit code {code}). See the console.')
        else:
            self.step_done.emit(label)
            self.index += 1
            self._next()

    def stop(self):
        if not self.proc:
            return
        self.stopped = True
        pid = int(self.proc.processId())
        self.line.emit('\n■ Stopping (process group %d)…' % pid)
        for sig, delay in ((signal.SIGTERM, 0), (signal.SIGKILL, 4000)):
            QTimer.singleShot(delay, lambda s=sig: self._signal(pid, s))

    def _signal(self, pid, sig):
        try:
            os.killpg(pid, sig)
        except (ProcessLookupError, PermissionError):
            pass
