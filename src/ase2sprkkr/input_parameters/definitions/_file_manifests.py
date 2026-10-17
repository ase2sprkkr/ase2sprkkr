"""Reusable pathname callbacks for declarative task files."""

import os


def _optional_value(section, name):
    try:
        return section[name]()
    except KeyError:
        return None


def tau_file(parameters):
    """Return the effective TAU-directory pathname for concrete parameters."""
    filename = _optional_value(parameters.TAU, "TAUDIR")
    if filename is None:
        filename = "{dataset}.tau.dir"
    else:
        filename = os.fspath(filename)
        if not filename:
            return filename
        if os.path.isabs(filename):
            return filename

    taupath = _optional_value(parameters.TAU, "TAUPATH")
    if taupath is None:
        return filename
    taupath = os.fspath(taupath)
    if not taupath:
        return taupath
    return os.path.join(taupath, filename)
