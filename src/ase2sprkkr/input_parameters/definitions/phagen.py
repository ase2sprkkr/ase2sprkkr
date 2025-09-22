""" PHAGEN task input parameters definition"""
# from ...common.grammar_types import
from .sections import CONTROL, TASK, TAU, ENERGY, SCF
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters \
    # ,InputValueDefinition as V

input_parameters = lambda: InputParameters(
      'phagen', [
        CONTROL('PHAGEN'),
        TASK('PHAGEN'),
        TAU,
        ENERGY(),
        SCF,
      ],
      info = "PHAGEN",
      description = "<TODO: change to something meaningfull>",
      executable = 'kkrscf',
      mpi=False
)
""" PHAGEN - PHAGEN task input parameters definition"""
