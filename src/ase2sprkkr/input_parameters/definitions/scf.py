""" SCF task input parameters definition"""
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V

input_parameters = InputParameters(
      'scf', [
        CONTROL('SCF').copy([
          V('KRWS', Integer(min=0, max=1), 1, info='If 0 RWS is taken from the potential file and scaled. If 1, RWS is calculated by scaling the muffin-tin radii by a common scaling factor. (This setting is forced in the case of FULLPOT.)'),
          V('KRMT', {
            '0' : 'RMT is taken from the potential file',
            '1' : 'RMT = min( x*RWS )',
            '2' : 'RMT = min( d_ij / 2 )',
            '3' : 'RMT from atomic charge density (=> KRWS=1)',
            '4' : 'RMT from atomic Hartree potential (=> KRWS=1)',
            '5' : 'RMT from total atomic potential (=> KRWS=1)',
            '6' : 'take average of 3 and 4 (=> KRWS=1)',
          }, default_value=None, info='It controls how the muffin-tin radii are calculated.',required=False)
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
