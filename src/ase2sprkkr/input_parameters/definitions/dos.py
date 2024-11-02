""" DOS task input parameters definition"""
from .sections import CONTROL, TAU, ENERGY, TASK, SITES
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters

input_parameters = lambda: InputParameters(
  'dos', [
      CONTROL('DOS'),
      TAU,
      ENERGY(
          emax = ( 1.0, 'value of the highest energy', None),
          defaults= {
            'GRID' : 3,
            'NE' : 300,
            'ImE' : 0.01,
          }),
      TASK('DOS'),
      SITES
  ],
  executable='kkrgen',
  mpi=True,
  info="DOS - The density of states computation"
)
""" DOS - density of states input parameters definition"""
