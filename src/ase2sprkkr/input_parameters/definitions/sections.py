""" Definitions of sections to be used in definitions of input parameters
(input files for SPR-KKR tasks) """

from ...common.grammar_types  import DefKeyword, SetOf, Flag, energy
from ..input_parameters_definitions import \
      SectionDefinition as Section, \
      ValueDefinition as V

def CONTROL(ADSI):
  """ Create the definition of the CONTROL section of the task input file.

  Parameters
  ----------
  ADSI: string
    the default value for the ADSI parameter of the resulting section


  Return
  ------
  CONTROL: SectionDefinition
  """
  return Section('CONTROL',[
      V('DATASET', str, 'case', required = True, help="Meaning of the parameter"),
      V('ADSI', DefKeyword(ADSI), required = True, help="Type of the computation -- do DFT selfconsistent cycle"),
      V('POTFIL', str, required=True, help="Potential file (see SPRKKR documentation for its format). It is not necessary to set it, it will be set by the calculator."),
      V('KRWS', int, required=False)
  ])

TAU = Section('TAU',[
      V('BZINT', DefKeyword('POINTS', 'ANOTHER_OPTION'), required=True),
      V('NKTAB', 250)
  ])
"""The definition of the TAU section of the task input file """

ENERGY = Section('ENERGY',[
      V('GRID', [5], required=True),
      V('NE', [32], required=True),
      V('ImE', energy, 0.0),
      V('EMIN', -0.2),
  ])
"""The definition of the ENERGY section of the task input file """

SCF = Section('SCF', [
      V('NITER', 200),
      V('MIX', 0.2),
      V('VXC', DefKeyword('VWN')),
      V('EFGUESS', 0.7),
      V('TOL', 0.00001),
      V('ISTBRY', 1)
  ])
"""The definition of the SCF section of the task input file """

SITES = Section('SITES', [
      V('NL', [3])
  ])
"""The definition of the SITES section of the task input file """

def TASK(TASK):
  """ Create the definition of the CONTROL section of the task input file.

  Parameters
  ----------
  TASK: string
    the default value for the TASK parameter of the resulting section

  Return
  ------
  TASK: SectionDefinition
  """
  return Section('TASK', [
    V('TASK', DefKeyword(TASK), name_in_grammar=False)
  ])

__all__ = [
    'CONTROL', 'TAU', 'ENERGY', 'SCF', 'SITES', 'TASK'
]
