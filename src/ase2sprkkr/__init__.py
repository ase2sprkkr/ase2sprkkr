"""
ase2sprkkr - ASE interface to SPR-KKR electron structure calculation.

This root package import a few most used objects to provide shortcuts to them.
"""

import warnings

__all__ = (
    "SPRKKR",
    "InputParameters",
    "Potential",
    "PotentialFileChanged",
    "FileState",
    "FileVersion",
    "FileVersionCaptureError",
    "SPRKKRAtoms",
    "promote_ase_atoms",
    "ResultFile",
    "TaskResult",
    "executables",
    "__version__",
    "config",
)

""" SPRKKR calculator to be used to calculate electronic structure using ASE """
from .sprkkr.calculator import SPRKKR  # NOQA: F401, E402

""" Input parameters object for SPRKKR calculation tasks """
from .input_parameters.input_parameters import InputParameters  # NOQA: F401, E402

""" Potential file object and its stable runtime exception contract """
from .potentials.potentials import Potential, PotentialFileChanged  # NOQA: F401, E402

""" Content-aware filesystem provenance shared with client applications """
from .common.file_versions import (  # NOQA: F401, E402
    FileState,
    FileVersion,
    FileVersionCaptureError,
)

""" An extension of ASE atoms object """
from .sprkkr.sprkkr_atoms import SPRKKRAtoms  # NOQA: F401, E402

""" For promoting ASE atoms to use ASE2SPRKKR extensions"""
promote_ase_atoms = SPRKKRAtoms.promote_ase_atoms

""" Structured SPR-KKR result file """
from .result_files import ResultFile  # NOQA: F401, E402

""" Read-only metadata for known SPR-KKR executables """
from .sprkkr.executables import executables  # NOQA: F401, E402


def __getattr__(name):
    if name == "OutputFile":
        warnings.warn(
            "ase2sprkkr.OutputFile is deprecated; use ase2sprkkr.ResultFile instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return ResultFile
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


""" Task results """
from .outputs.task_result import TaskResult  # NOQA E402

""" Version of the package """
from .version import __version__  # NOQA: F401, E402


def _init():
    import ase  # NOQA: F401
    from .ase.register import register
    from .configuration import load_user_preferences  # noqa: F401

    register()
    import os

    if not os.environ.get("ASE2SPRKKR_NO_USER_PROFILE", False):
        load_user_preferences()


_init()
del _init

from .configuration import config as config
