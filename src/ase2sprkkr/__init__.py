"""
ase2sprkkr - ASE interface to SPR-KKR electron structure calculation.

This root package import a few most used objects to provide shortcuts to them.
"""

""" SPRKKR calculator to be used to calculate electronic structure using ASE """
from .sprkkr.calculator import SPRKKR

""" Input parameters object for SPRKKR calculation tasks """
from .input_parameters.input_parameters import InputParameters

""" Potential file object for SPRKKR calcualtions """
from .potentials.potentials import Potential

""" An extension of ASE atoms object """
from .sprkkr.sprkkr_atoms import SPRKKRAtoms

""" SPRKKR Output File """
from .outputs.output_files import OutputFile

""" Version of the package """
from .version import __version__

def _init():
   import ase
   from ase.calculators.calculator import register_calculator_class
   register_calculator_class('sprkkr', SPRKKR)

_init()
del _init
