""" DOS task input parameters definition"""
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V
from ...common.doc import process_input_parameters_definition

input_parameters = InputParameters(
      'jxc', [
          CONTROL('JXC'),
          TAU,
          MODE,
          STRCONST,
          ENERGY,
          TASK('JXC').copy([
              V('CLURAD', 2.2, info=
    """The radius of a sphere restricting the cluster
    of atoms for which the exchange coupling
    constants will be calculated"""
                )
          ]),
          SITES
      ],
    executable='kkrgen',
    mpi=True
)
""" DOS task input parameters definition"""

process_input_parameters_definition(__name__)
