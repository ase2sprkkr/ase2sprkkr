#!/bin/bash
set -e

PYTHON_MINOR=$(python3 -c 'import sys; print(sys.version_info.minor)')

if (( PYTHON_MINOR >= 14 )); then
    echo "Python >= 3.14 detected. Installing monty 2025.3.3 ignoring Python version metadata..."
    pip install --ignore-requires-python monty==2025.3.3 pymatgen==2025.10.7
    # The fork carries Python 3.14 compatibility changes.  Do not install it
    # on older Python versions: its build requires NumPy 2.x, which is not
    # available for Python 3.8 (and is unnecessary there).
    python3 -m pip install git+https://github.com/lokik/spglib.git
fi

python3 -m pip install .
