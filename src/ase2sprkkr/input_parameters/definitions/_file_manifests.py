"""Helpers for task file manifests."""

import os


def resolve_file(directory, filename):
    """Resolve a file name relative to the calculation directory."""
    filename = os.fspath(filename)
    return filename if os.path.isabs(filename) else os.path.join(directory, filename)


def _optional_value(section, name):
    try:
        return section[name]()
    except KeyError:
        return None


def tau_file(parameters, directory):
    """Return the effective native TAU-directory file path."""
    filename = _optional_value(parameters.TAU, "TAUDIR")
    if not filename:
        filename = f"{parameters.CONTROL.DATASET()}.tau.dir"

    filename = os.fspath(filename)
    if os.path.isabs(filename):
        return filename

    taupath = _optional_value(parameters.TAU, "TAUPATH")
    base = resolve_file(directory, taupath) if taupath else directory
    return os.path.join(base, filename)


def bsf_file(parameters, directory):
    """Return the native BSF output path for the concrete dataset/ADSI."""
    dataset = parameters.CONTROL.DATASET()
    adsi = parameters.CONTROL.ADSI()
    if not dataset or not adsi:
        return None
    return resolve_file(directory, f"{dataset}_{adsi}.bsf")
