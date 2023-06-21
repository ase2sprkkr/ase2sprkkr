""" SCF task input parameters definition"""
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V

input_parameters = InputParameters(
      'scf', [
        CONTROL('SCF').copy([
          V('KRWS', Integer(min=0, max=1), 1),
          V('KRMT', Integer(min=0, max=6), required=False)
        ]),
        TAU,
        ENERGY,
        SCF,
        SITES,
        STRCONST,
        CPA,
        MODE
      ],
      info = "SCF - calculate a self-consistent potential",
      description = "",
      executable = 'kkrscf',
      mpi = True
)
""" SCF task input parameters definition"""

from ...common.doc import process_input_parameters_definition
process_input_parameters_definition(__name__)
