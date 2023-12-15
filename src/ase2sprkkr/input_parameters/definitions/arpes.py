""" ARPES task input parameters definition """
from ...common.grammar_types import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V, \
    InputSectionDefinition as Section
from ...sprkkr.sprkkr_grammar_types import Site

input_parameters = InputParameters(
    'arpes', [
    CONTROL('ARPES'),
    TAU,
    ENERGY.copy([
        V('EMINEV', -10., info='Minimum of the energy window in eV with respect to the Fermi level'),
        V('EMAXEV', -10., info='Maximum of the energy window in eV with respect to the Fermi level'),
        V('EWORK_EV', 4.2, info='Inner potential of the bulk crystal in eV'),
        V('IMV_INI_EV', 0.05, info='Imaginary part of the potential in eV (initial state)'),  # alternatively you can use VIL (in eV) or IMV_INI (in Ry)'),
        V('IMV_FIN_EV', 5., info='Imaginary part of the potential in eV (final state)'),  # alternatively you can use VIH (in eV) or IMV_FIN (in Ry)'),
      ],
      remove = ['EMIN'],
      defaults = { 'GRID' : 1, 'NE' : 300 }
    ),

    SITES.copy(defaults = {'NL' : 4 }),

    Section('TASK', [
      V('TASK', DefKeyword({
        'ARPES' : 'Angle resolved photoemission spectroscopy',
        'AIPES' : 'Angle integrated photoemission spectroscopy',
        'SPLEED' : 'Spin polarized LEED (experimental feature)',
        'BAND'  : 'band structure calculations (experimental feature)'
      }), name_in_grammar=False, info='Type of the calculation'),
      V('IQ_AT_SURF', Site.I, 1),
      V('MILLER_HKL', SetOf(int, length=3), [0,0,1]),
      V('CRYS_VEC', True, info='Miller indices with respect to crystalographic primitive vectors'),

      V('STRVER', 1, is_expert=True, is_always_added=True, info="Set to 0 to supply the ARPES input file 'struc.inp' manually (and do not generate it)."),
      V('INPVER', 1, is_expert=True, is_always_added=True, info='Set to 0 to use an old input.inp from old rslab'),
    ]),

    Section('SPEC_PH', [
      V('THETA', 45., info='Direction of the photon (the polar coordinate)'),
      V('PHI', 0., info='Direction of the photon (the azimuth coordinate)'),
      V('POL_P', DefKeyword('P', 'S', 'C+', 'C-'), info='Polarization of the light'),
      V('EPHOT', 25., info='Photon energy in eV'),
      #Expert
      V('ALQ', expert=45., info='Alignment of polarization vector or pol.ellipsis'),
      V('DELQ', expert=0., info='Phase shift between real and imaginary part of e-vector, delq=90 for circular polarized light'),
      V('NPOL', Keyword({
        0: 'unpolarized and p-s dichroism for the calculation',
        1: 'p-pol or rcp or elliptical (depends on icirc, etc.)',
        2: 's-pol or lcp or elliptical (depends on icirc, etc.)',
        3: 'dichroism (ddad, ldad)'
        }), expert=1, info='Controls the polarization and dichroism'),
      V('ICIRC', Keyword({
        0: 'elliptically pol. light: alq, delq arbitrary',
        1: 'linear pol. light: alq arbitrary, delq=0',
        2: 'circular pol. light: alq=45, delq = 90'
      }), expert=1, info='controls the polarization and dichroism'),
      V('IDREH', Keyword({
          0: 'linearly polarized (equals icirc=1)',
          1: 'right circular polarization',
         -1: 'left circular polarization'
         }), expert=0, info='Helicity of the photons'),
      V('IFSP', Keyword({
          0: 'fixed',
          1: 'variable'
          }), is_expert=True, required=False, info='Photon azimuth angle type'),
      V('THETA_FIX', float, is_expert=True, info='Light and electrons are at fixed polarization angle')
    ], info=''),

    Section('SPEC_EL', [
      V('THETA', Range(float), 45., info='(A range of) scattering angle(s) (the azimuth coordinate)'),
      V('PHI', Range(float), 0., info='(A range of) scattering angle(s) (the polar coordinate)'),
      V('NT', 1, info='Number of angular values for a rotation in the polar coordinate.'),
      V('NP', 1, info='Number of angular values for a rotation in the azimuth coordinate.'),
      V('POL_E', DefKeyword('PZ')),
      V('SPOL', int, required=False),
      V('PSPIN', SetOf(float, length=3), required=False),
      V('BETA1', float, required=False, info='Begin of the rotation'),
      V('BETA2', float, required=False, info='End of the rotation'),
      V('ROTAXIS',SetOf(int, length=3),required=False, info='Axis of the rotation'),
      #expert
      V('TYP', Keyword({0: "i(e) diagram",
                        1: "rotation diagram -> phi scan",
                        2: "scattering-angle diagram -> theta scan",
                        3: "orthonormal projection",
                        4: "stereographic projection"},
                       description = '3,4 only for angular resolved\npe (ups, xps) note: nt=np-> nx,ny'
              ), expert=1,
              info='Crystal coordinats in splout, xpsrun, or upsrun'),
      V('ISTR', Array(int, length=2), expert=[0,0], info="beam number (h,k)"),
      V('POL0', Array(int, length=3), expert=[0,0,0], info="initial pol."),
      V('POL0L', Array(int, length=3), expert=[0,0,0], info="initial pol. in the laboratory system"),
      V('Q1', complex, expert=1. + 0.j, info="Amplitude 1 of the photoelectron used in spin polarized calculations"),
      V('Q2', complex, expert=0. + 0.j, info="Amplitude 2 of the photoelectron used in spin polarized calculations"),
      V('Q3', complex, expert=0. + 0.j, info="Amplitude 3 of the photoelectron used in spin polarized calculations"),
      V('Q4', complex, expert=1. + 0.j, info="Amplitude 4 of the photoelectron used in spin polarized calculations"),
    ]),

    Section('SPEC_STR', [
      V('N_LAYDBL', SetOf(int), [10,10]),
      V('NLAT_G_VEC', 57),
      V('N_LAYER', 50),
      V('SURF_BAR', SetOf(float), [0.25,0.25]),
      V('TRANSP_BAR', flag, False)
    ]),
  ],
  executable = 'kkrspec',
  info = 'ARPES - angle resolved photoemission spectroscopy'
)
""" ARPES task input parameters definition"""

from ...common.doc import process_input_parameters_definition
process_input_parameters_definition(__name__)
