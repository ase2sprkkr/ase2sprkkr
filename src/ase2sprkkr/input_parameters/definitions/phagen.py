""" PHAGEN task input parameters definition"""
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V

input_parameters = InputParameters(
      'phagen', [
        CONTROL('PHAGEN'),
        TAU,
        ENERGY,
        SCF,
      ],
      info = "PHAGEN",
      description = "<TODO: change to something meaningfull>",
      executable = 'kkrscf',
      mpi=False
)
""" PHAGEN - PHAGEN task input parameters definition"""

from ...common.doc import process_input_parameters_definition
process_input_parameters_definition(__name__)
