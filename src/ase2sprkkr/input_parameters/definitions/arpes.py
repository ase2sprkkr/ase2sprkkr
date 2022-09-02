""" ARPES task input parameters definition """
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    ValueDefinition as V, \
    SectionDefinition as Section

#: ARPES task input parameters definition
input_parameters = InputParameters(
    'arpes', [
    CONTROL('ARPES'),
    TAU,
    ENERGY.copy([
        V('EMINEV', -10., help='Minimum of the energy window in eV with respect to the Fermi level'),
        V('EMAXEV', -10.),
        V('EWORKEV', 4.2),
        V('IMV_INI_EV', 0.05),
        V('IMV_FIN_EV', 5.),
      ],
      remove = ['EMIN'],
      defaults = { 'GRID' : 1, 'NE' : 300 }
    ),
    SITES.copy(defaults = {'NL' : 4 }),
    TASK('ARPES').copy([
      V('STRVER', 0),
      V('IQ_AT_SURF', 2),
      V('MILLER_HKL', SetOf(int, length=3), [0,0,1]),
      V('NTMP', 1),
      V('TMPMIN', 11.),
      V('CTMPMAX', 11.),
      V('CTMPMAX', 11.),
      V('VIBRA', flag, True),
      V('CNVIBRA', 14),
    ]),
    Section('SPEC_PH', [
      V('THETA', 45.),
      V('PHI', 0.),
      V('POL_P', 'P'),
      V('EPHOT', 6675.),
    ]),
    Section('SPEC_EL', [
      V('THETA', Range(float), 45.),
      V('PHI', Range(float), 0),
      V('NT', 1),
      V('NP', 1),
      V('POL_E', DefKeyword('PZ')),
      V('SPOL', int, required=False),
      V('PSPIN', SetOf(float, length=3), required=False)
    ]),
    Section('SPEC_STR', [
      V('N_LAYDBL', SetOf(int), [10,10]),
      V('NLAT_G_VEC', 57),
      V('N_LAYER', 120),
      V('SURF_BAR', SetOf(float), [0.25,0.25]),
      V('TRANSP_BAR', flag, False)
    ]),
  ],
  executable = 'kkrspec'
)
""" ARPES task input parameters definition"""
