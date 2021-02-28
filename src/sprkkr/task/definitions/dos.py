from .option_types  import *
from .configuration_definitions import \
    TaskDefinition as Task, \
    SectionDefinition as Section, \
    OptionDefinition as O

task = Task(
      'dos', [
      CONTROL('DOS'),
      TAU,
      ENERGY.copy([
          O('EMAX', 1)
      ], defaults= {
          'GRID' : 5
      }),
      TASK('DOS'),
      SITES
])
