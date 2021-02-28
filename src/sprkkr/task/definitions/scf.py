from .option_types  import *
from .sections import *
from .configuration_definitions import \
    TaskDefinition as Task, \
    OptionDefinition as O

task = Task(
      'scf', [
        CONTROL('SCF').copy([
          O('KRWS', 1)
        ]),
        TAU,
        ENERGY.copy( defaults = {'NE' : 32}),
        SCF,
        SITES
      ],
      description = "SCF calculation",
      help = "Set at least options ..., ...., ....."
)
