""" DOS task input parameters definition"""
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    ValueDefinition as V \

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
""" DOS task input parameters definition"""
