def _sections():
  from ...common.grammar_types  import Keyword, DefKeyword, SetOf, Flag, energy
  from ..task_definitions import \
      SectionDefinition as Section, \
      ValueDefinition as V

  sections = {}
  sections['CONTROL'] = lambda x: Section('CONTROL',[
      V('DATASET', str, 'case', required = True, help="Meaning of the parameter"),
      V('ADSI', DefKeyword(x), required = True, help="Type of the computation -- do DFT selfconsistent cycle"),
      V('POTFIL', str, required=True, help="Potential file (see SPRKKR documentation for its format). It is not necessary to set it, it will be set by the calculator."),
      V('KRWS', int)
  ])

  sections['TAU'] = Section('TAU',[
      V('BZINT', DefKeyword('POINTS', 'ANOTHER_OPTION'), required=True),
      V('NKTAB', 250)
  ])

  sections['ENERGY'] = Section('ENERGY',[
      V('GRID', [5], required=True),
      V('NE', [300], required=True),
      V('ImE', energy, 0.0),
      V('EMIN', -0.2),
  ])

  sections['SCF'] = Section('SCF', [
      V('NITER', 200),
      V('MIX', 0.2),
      V('VXC', DefKeyword('VWN')),
      V('EFUESS', 0.7),
      V('TOL', 0.00001),
      V('ISTBRY', 1)
  ])

  sections['SITES'] = Section('SITES', [
      V('NL', [3])
  ])

  sections['TASK'] = lambda x: Section('TASK', [
    V('TASK', DefKeyword(x))
  ])

  return sections

locals().update(_sections())
del _sections