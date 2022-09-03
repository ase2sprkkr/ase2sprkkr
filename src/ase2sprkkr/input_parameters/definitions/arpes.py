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
        V('EMAXEV', -10., help='Maximum of the energy window in eV with respect to the Fermi level'),
        V('EWORKEV', 4.2, help='Inner potential of the bulk crystal in eV'),
        V('IMV_INI_EV', 0.05, help='Imaginary part of the potential in eV (initial state)'), # alternatively you can use VIL (in eV) or IMV_INI (in Ry)'),
        V('IMV_FIN_EV', 5., help='Imaginary part of the potential in eV (final state)'), # alternatively you can use VIH (in eV) or IMV_FIN (in Ry)'),
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
      V('POL_P', DefKeyword('P', 'S', 'C+', 'C-'), help='Polarization of the light'),
      V('EPHOT', 6675., help='Photon energy in eV'),
      #Expert
      V('ALQ', expert=45., help='Alignment of polarization vector or pol.ellipsis'),
      V('DELQ', expert=0., help='Phase shift between real and imaginary part of e-vector, delq=90 for circular polarized light'),
      V('NPOL', Keyword({
        0: 'unpolarized and p-s dichroism for the calculation',
        1: 'p-pol or rcp or elliptical (depends on icirc, etc.)',
        2: 's-pol or lcp or elliptical (depends on icirc, etc.)',
        3: 'dichroism (ddad, ldad)'
        }), expert=1, help='Controls the polarization and dichroism'),
      V('ICIRC', Keyword({
        0: 'elliptically pol. light: alq, delq arbitrary',
        1: 'linear pol. light: alq arbitrary, delq=0',
        2: 'circular pol. light: alq=45, delq = 90'
      }), expert=1, help='controls the polarization and dichroism'),
      V('IDREH', Keyword({
          0: 'linearly polarized (equals icirc=1)',
          1: 'right circular polarization',
         -1: 'left circular polarization'
         }), expert=0, help='Helicity of the photons'),
      V('IFSP', Keyword({
          0: 'fixed',
          1: 'variable'
          }), is_expert=True, required=False, help='Photon azimuth angle type'),
      V('THETA_FIX', float, is_expert=True, help='Light and electrons are at fixed polarization angle')
    ]),

    Section('SPEC_EL', [
      V('THETA', Range(float), 45., help='Scattering angle'),
      V('PHI', Range(float), 0, help='Scattering angle'),
      V('NT', 1, help='Number of angular values for a rotation in polar coordinate.'),
      V('NP', 1, help='Number of angular values for a rotation in azimuth coordinate.'),
      V('POL_E', DefKeyword('PZ')),
      V('SPOL', int, required=False),
      V('PSPIN', SetOf(float, length=3), required=False),
      #expert
      V('TYP', Keyword({0: "i(e) diagram",
                        1: "rotation diagram -> phi scan",
                        2: "scattering-angle diagram -> theta scan",
                        3: "orthonormal projection",
                        4: "stereographic projection"},
                      additional_description = '3,4 only for angular resolved\npe (ups, xps) note: nt=np-> nx,ny'
              ), expert=1,
              help='Crystal coordinats in splout, xpsrun, or upsrun'),
      V('ISTR', Array(int, length=2), expert=[0,0], help="beam number (h,k)"),
      V('POL0', Array(int, length=3), expert=[0,0,0], help="initial pol."),
      V('POL0L', Array(int, length=3), expert=[0,0,0], help="initial pol. in the laboratory system"),
      V('Q1', complex, expert=1.+0.j, help="Amplitude 1 of the photoelectron used in spin polarized calculations"),
      V('Q2', complex, expert=0.+0.j, help="Amplitude 2 of the photoelectron used in spin polarized calculations"),
      V('Q3', complex, expert=0.+0.j, help="Amplitude 3 of the photoelectron used in spin polarized calculations"),
      V('Q4', complex, expert=1.+0.j, help="Amplitude 4 of the photoelectron used in spin polarized calculations"),
    ]),

    Section('SPEC_STR', [
      V('N_LAYDBL', SetOf(int), [10,10]),
      V('NLAT_G_VEC', 57),
      V('N_LAYER', 120),
      V('SURF_BAR', SetOf(float), [0.25,0.25]),
      V('TRANSP_BAR', flag, False)
    ]),
  ],
  executable = 'kkrspec',
  help = 'ARPES - angle resolved photoemission spectroscopy'
)
""" ARPES task input parameters definition"""
