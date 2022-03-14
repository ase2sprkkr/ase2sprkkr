""" PHAGEN task input parameters definition"""
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    ValueDefinition as V

input_parameters = InputParameters(
      'phagen', [
        CONTROL('PHAGEN'),
        TAU,
        ENERGY,
        SCF,
      ],
      description = "SCF calculation",
      help = "Set at least options ..., ...., .....",
      executable = 'kkrscf',
      mpi=False
)
""" PHAGEN task input parameters definition"""
