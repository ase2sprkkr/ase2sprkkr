from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    ValueDefinition as V

input_parameters = InputParameters(
      'phagen', [
        CONTROL('PHAGEN'),
        TAU,
        ENERGY.copy( defaults = {'NE' : 32}),
        SCF,
      ],
      description = "SCF calculation",
      help = "Set at least options ..., ...., .....",
      executable = 'kkrgen'
)
