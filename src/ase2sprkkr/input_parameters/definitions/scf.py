""" SCF task input parameters definition"""
from .sections import CONTROL, TAU, ENERGY, SCF, SITES, STRCONST, CPA, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters \
    #   ,InputValueDefinition as V
from ...common.doc import process_input_parameters_definition

input_parameters = InputParameters(
      'scf', [
        CONTROL('SCF'),
        TAU,
        ENERGY,
        SCF,
        SITES,
        STRCONST,
        CPA,
        MODE
      ],
      info = "SCF - calculate a self-consistent potential",
      description = "",
      executable = 'kkrscf',
      mpi = True
)
""" SCF task input parameters definition"""

process_input_parameters_definition(__name__)
