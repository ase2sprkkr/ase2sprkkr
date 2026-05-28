"""SCF task input parameters definition"""

from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE, TASK
from ..input_parameters_definitions import InputParametersDefinition as InputParameters
#   ,InputValueDefinition as V

def input_parameters(): return InputParameters(
    "scf",
    [CONTROL("SCF"), TASK("SCF"), TAU, ENERGY(), SCF, SITES, STRCONST, CPA, MODE],
    info="SCF - Calculate a self-consistent potential",
    description="",
    executable="kkrscf",
    mpi=True,
)
""" SCF task input parameters definition"""
