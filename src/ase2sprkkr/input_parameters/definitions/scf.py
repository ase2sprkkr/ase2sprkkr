from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    ValueDefinition as V

input_parameters = InputParameters(
      'scf', [
        CONTROL('SCF').copy([
          V('KRWS', 1)
        ]),
        TAU,
        ENERGY.copy( defaults = {'NE' : 32}),
        SCF,
        SITES
      ],
      description = "SCF calculation",
      help = "Set at least options ..., ...., .....",
      executable = 'kkrscf',
      mpi = True
)