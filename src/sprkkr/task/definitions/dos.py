from ...common.grammar_types  import *
from .sections import *
from ..task_definitions import \
    TaskDefinition as Task, \
    ValueDefinition as V \

task = Task(
      'dos', [
      CONTROL('DOS'),
      TAU,
      ENERGY.copy([
          V('EMAX', 1.0)
      ], defaults= {
          'GRID' : 5
      }),
      TASK('DOS'),
      SITES
])
