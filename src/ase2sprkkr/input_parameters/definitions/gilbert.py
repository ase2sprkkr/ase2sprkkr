""" Gilbert task input parameters definition """
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V, \
    InputSectionDefinition as Section

input_parameters = InputParameters(
    'gilbert', [
    CONTROL('GILBERT'),
    TAU,
    ENERGY.copy(remove = ['EMIN'],
      defaults = { 'GRID' : 3, 'NE' : 1 }
    ),
    SITES.copy(defaults = {'NL' : 4 }),

    TASK('Gilbert').copy([
      V('NTMP', int, required=False, info='Number of temperature points used for α(T)'),
      V('SETFLUCT', Keyword({
        'MLIN': 'Use linear temperature grid. Takes account only the electron scattering due to lattice vibrations.',
        'M_T': 'If NFTET = 1 and NFPHI = 1 but NVIBRA is bigger than 1, '
               'only lattice vibrations are taken into account. If NVIBRA = 1 but '
               'NFTET is bigger than 1 and NFPHI is bigger than 1, only spin '
               'fluctuations are taken into account.',
      }), default_value=None, required=False, info='Finite temperature calculation mode'),
      V('NVIBRA', int, required=False, info='Number of directions for atomic displacements'
                                            ' representing thermal lattice vibrations'),
      V('TMPMIN', real, required=False, info='Lower limit of the temperature region (SETLFUNC=MLIN).'),
      V('TMPMAX', real, required=False, info='Upper limit of the temperature region (SETLFUNC=MLIN)'),
      V('FLUCTFIL', str, required=False, info='Data ﬁle which contains the information about '
                                              'temperature dependent magnetizstion, taken from the '
                                              'experiment or Monte Carlo simulations'),
      V('NFTET', int, required=False, info='Number of grid points specifying θ angle (SETLFUNC=M_T)'),
      V('NFPHI', int, required=False, info='Number of grid points specifying φ angle (SETLFUNC=M_T)'),
    ]),

 ],
  executable = 'kkrchi',
  info = 'GILBERT - The Gilbert damping parameter calculation'
)
""" The Gilbert damping input parameters definition"""

from ...common.doc import process_input_parameters_definition
process_input_parameters_definition(__name__)
