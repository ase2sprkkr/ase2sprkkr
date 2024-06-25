"""
ase2sprkkr - ASE interface to SPR-KKR electron structure calculation.

This root package import a few most used objects to provide shortcuts to them.
"""

""" Configuration """
from . import config

""" SPRKKR calculator to be used to calculate electronic structure using ASE """
from .sprkkr.calculator import SPRKKR  # NOQA: E402

""" Input parameters object for SPRKKR calculation tasks """
from .input_parameters.input_parameters import InputParameters  # NOQA: F401, E402

""" Potential file object for SPRKKR calcualtions """
from .potentials.potentials import Potential  # NOQA: F401, E402

""" An extension of ASE atoms object """
from .sprkkr.sprkkr_atoms import SPRKKRAtoms  # NOQA: F401, E402

""" SPRKKR Output File """
from .output_files.output_files import OutputFile  # NOQA: F401, E402

""" Version of the package """
from .version import __version__  # NOQA: F401, E402


def _init():
    import ase   # NOQA: F401
    import platformdirs
    import os

    from ase.calculators.calculator import register_calculator_class
    register_calculator_class('sprkkr', SPRKKR)

    file = os.path.join(platformdirs.user_config_dir('ase2sprkkr', 'ase2sprkkr'), '__init__.py')
    try:
       if os.path.isfile(file):
           import types
           import importlib.machinery
           loader = importlib.machinery.SourceFileLoader('ase2sprkkr.personal', file)
           mod = types.ModuleType(loader.name)
           loader.exec_module(mod)
    except Exception as e:
        import warnings
        warnings.warn(f'Can not import {file} file with the user preferences: \n{e}')


_init()
del _init
