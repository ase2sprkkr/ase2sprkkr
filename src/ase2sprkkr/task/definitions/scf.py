from ...common.grammar_types  import *
from .sections import *
from ..task_definitions import \
    TaskDefinition as Task, \
    ValueDefinition as V

task = Task(
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
