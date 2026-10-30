#!/usr/bin/env bash
# FDTD_meep v1.2.0 | 2026-10-30 | Carlos Chacón
# Installer (macOS / Linux): creates ./.env with Meep + MPI + PyQt5 and runs the tests.
set -euo pipefail
cd "$(dirname "$0")"

if command -v micromamba >/dev/null 2>&1; then
    MM=micromamba
elif command -v mamba >/dev/null 2>&1; then
    MM=mamba
elif command -v conda >/dev/null 2>&1; then
    MM=conda
elif [[ "$(uname)" == "Darwin" ]] && command -v brew >/dev/null 2>&1; then
    echo "Installing micromamba with Homebrew..."
    brew install micromamba
    MM=micromamba
else
    echo "micromamba, mamba or conda is required (Meep is distributed through conda-forge)."
    echo "Install micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
    exit 1
fi

if [[ -x ./.env/bin/python ]]; then
    echo "Updating the existing environment ./.env with $MM..."
    if [[ "$MM" == micromamba ]]; then
        micromamba install -y -p ./.env -f environment.yml
    else
        "$MM" env update -p ./.env -f environment.yml
    fi
else
    echo "Creating the environment ./.env with $MM (Meep, MPI, PyQt5; takes a few minutes)..."
    if [[ "$MM" == micromamba ]]; then
        micromamba create -y -p ./.env -f environment.yml
    else
        "$MM" env create -p ./.env -f environment.yml
    fi
fi
./.env/bin/python -m pip install --quiet "miepython==3.3.0"
chmod +x run_gui.sh

echo "Running the tests..."
./.env/bin/python -m unittest discover -s tests

echo
echo "Done. Start the interface with:   ./run_gui.sh"
