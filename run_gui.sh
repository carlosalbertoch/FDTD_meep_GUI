#!/usr/bin/env bash
# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Starts the graphical interface:  ./run_gui.sh [--config inputs/config/<file>.json]
cd "$(dirname "$0")"
if [[ ! -x ./.env/bin/python ]]; then
    echo "Environment not found. Run ./install.sh first."
    exit 1
fi
exec ./.env/bin/python -m fdtd_meep_gui "$@"
