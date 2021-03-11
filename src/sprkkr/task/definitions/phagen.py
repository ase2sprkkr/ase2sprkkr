from ...common.grammar_types  import *
from .sections import *
from ..task_definitions import \
    TaskDefinition as Task, \
    ValueDefinition as V

task = Task(
      'phagen', [
        CONTROL('PHAGEN'),
        TAU,
        ENERGY.copy( defaults = {'NE' : 32}),
        SCF,
      ],
      description = "SCF calculation",
      help = "Set at least options ..., ...., ....."
)
