"""
ase2sprkkr - ASE interface to SPR-KKR electron structure calculation.

This root package import a few most used objects to provide shortcuts to them.
"""

import warnings

""" SPRKKR calculator to be used to calculate electronic structure using ASE """
from .sprkkr.calculator import SPRKKR  # NOQA: F401, E402

""" Input parameters object for SPRKKR calculation tasks """
from .input_parameters.input_parameters import InputParameters  # NOQA: F401, E402

""" Potential file object for SPRKKR calcualtions """
from .potentials.potentials import Potential  # NOQA: F401, E402

""" An extension of ASE atoms object """
from .sprkkr.sprkkr_atoms import SPRKKRAtoms  # NOQA: F401, E402

""" For promoting ASE atoms to use ASE2SPRKKR extensions"""
promote_ase_atoms = SPRKKRAtoms.promote_ase_atoms

""" Structured SPR-KKR result file """
from .result_files import ResultFile  # NOQA: F401, E402


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

""" Enhance the atoms with sites array and other ASE2SPRKKR enhancements """
promote_ase_atoms = SPRKKRAtoms.promote_ase_atoms

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
