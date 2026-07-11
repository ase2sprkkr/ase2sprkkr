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

# Keep the wheel long enough to verify that CMake actually installed both
# compiled extensions into it.  ``pip install .`` hides its ephemeral wheel,
# so a successful build alone does not prove that the extensions were packed.
WHEEL_DIR=$(mktemp -d)
python3 -m pip wheel --no-deps --wheel-dir "$WHEEL_DIR" .
python3 - "$WHEEL_DIR" <<'PY'
import pathlib
import sys
import zipfile

wheel_dir = pathlib.Path(sys.argv[1])
wheels = list(wheel_dir.glob("ase2sprkkr-*.whl"))
if len(wheels) != 1:
    raise RuntimeError(f"Expected exactly one ASE2SPRKKR wheel, found: {wheels}")

with zipfile.ZipFile(wheels[0]) as wheel:
    names = wheel.namelist()
    for module in ("spheres", "symmetry"):
        prefix = f"ase2sprkkr/bindings/xband/{module}."
        matches = [
            name for name in names
            if name.startswith(prefix) and name.endswith((".so", ".pyd"))
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"Wheel {wheels[0]} does not contain exactly one compiled "
                f"{module} extension; found: {matches}"
            )
        print(f"Verified wheel extension: {matches[0]}")
PY

python3 -m pip install "$WHEEL_DIR"/ase2sprkkr-*.whl
