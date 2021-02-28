def _sections():
  from .option_types  import Integer, Unsigned, Keyword, SetOf, Flag
  from .configuration_definitions import \
      TaskDefinition as Task, \
      SectionDefinition as Section, \
      OptionDefinition as O

  sections['CONTROL'] = lambda x: Section('CONTROL',[
      O('ADSI', Keyword(x), required = True, help="Type of the computation -- do DFT selfconsistent cycle"),
      O('DATASET', str, required = True, help="Meaning of the parameter"),
      O('POTFIL', str, required=True, help="Potential file (see SPRKKR documentation for its format)"),
      O('KRWS', int)
  ])

  sections['TAU'] = Section('TAU',[
      O('BZINT', Keyword('POINTS', 'ANOTHER_OPTION'), required=True),
      O('NKTAB', 250)
  ])

  sections['ENERGY'] = Section('ENERGY',[
      O('GRID', [5], required=True),
      O('NE', [300], required=True),
      O('ImE', 0.0),
      O('EMIN', -0.2),
  ])

  sections['SCF'] = Section('SCF', [
      O('NITER', 200),
      O('MIX', 0.2),
      O('VXC', Keyword('VWN')),
      O('EFUESS', 0.7),
      O('TOL', 0.00001),
      O('ISTBRY', 1)
  ])

  sections['SITES'] = Section('SITES', [
      O('NL', [3])
  ])

  sections['TASK'] = lambda x: Section('TASK', [
    O('TASK', Keyword(x))
  ])

  return sections

locals.update(_sections())
del _sections
