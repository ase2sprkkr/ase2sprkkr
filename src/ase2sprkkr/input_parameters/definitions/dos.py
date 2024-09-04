""" DOS task input parameters definition"""
from .sections import CONTROL, TAU, ENERGY, TASK, SITES
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V
from ...common.doc import process_input_parameters_definition

input_parameters = InputParameters(
  'dos', [
      CONTROL('DOS'),
      TAU,
      ENERGY.copy([
          V('EMAX', 1.0)
      ], defaults= {
          'GRID' : 3,
          'NE' : 300,
          'ImE' : 0.01,
      }),
      TASK('DOS'),
      SITES
  ],
  executable='kkrgen',
  mpi=True
)
""" DOS - density of states input parameters definition"""

process_input_parameters_definition(__name__)
