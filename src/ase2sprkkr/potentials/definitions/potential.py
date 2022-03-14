""" Here, the format of the potential_file is defined """

from .sections import *

def fce():
  from ...common.misc import OrderedDict
  from  ...common.grammar_types import \
        Keyword, DefKeyword, Table, Array, Sequence, \
        integer, Date, boolean, line_string
  from ..potential_definitions import \
        PotSectionDefinition, \
        ASEArraySectionDefinition, \
        PotValueDefinition as V,\
        PotentialDefinition, \
        Separator
  import datetime

  sections = []
  def Section(*args, cls=PotSectionDefinition, **kwargs):
     x=cls(*args, **kwargs)
     sections.append(x)
     return x

  def ArraySection(array_name):
      def factory(*args, **kwargs):
          return ASEArraySectionDefinition(*args, array_name = array_name, **kwargs)
      return factory

  Section('HEADER', [
    Separator(),
    V('HEADER', line_string, lambda x: f'SPR-KKR potential file, created at {datetime.datetime.now()}'),
    Separator(),
    V('TITLE', line_string, 'Created by ASE-SPR-KKR wrapper'),
    V('SYSTEM', line_string, lambda x: 'System: {}'.format(x._get_root_container().atoms.symbols) ),
    V('PACKAGE', line_string, 'SPR-KKR'),
    V('FORMAT', Sequence(int, Date(prefix='(', postfix=')'), names = ['VERSION', 'DATE']),
          default_value = [7, datetime.datetime(2007,5,21)]),
  ], name_in_grammar = False)

  Section('GLOBAL SYSTEM PARAMETER', cls = GlobalSystemParameterDefinition)

  Section('SCF-INFO', [
    V('INFO', line_string, 'NONE'),
    V('SCFSTATUS', DefKeyword('START', 'CONVERGED', 'ITR-BULK')),
    V('FULLPOT', False),
    V('BREITINT', False),
    V('NONMAG', False, alternative_names='NOMAG'),
    V('ORBPOL', str, 'NONE'),
    V('EXTFIELD', False),
    V('BLCOUPL', False),
    V('BEXT', 0.0),
    V('SEMICORE', False),
    V('LLOYD', False),
    V('NE', Array(int), is_optional = True),
    V('IBZINT', int, is_optional = True),
    V('NKTAB', int, is_optional = True),
#   V('XC-POT', str, is_optional = True),
    V('SCF-ALG', Keyword('BROYDEN', 'BROYDEN2'), is_optional = True),
    V('SCF-ITER', 0),
    V('SCF-MIX', 0.2),
    V('SCF-TOL', 1e-5),
    V('RMSAVV', 999999.),
    V('RMSAVB', 999999.),
    V('EF', 999999.),
    V('VMTZ', 0.7),
  ], )

  Section('LATTICE', cls = LatticeSectionDefinition)
  Section('SITES',   cls = SitesSectionDefinition)
  Section('OCCUPATION',   cls = OccupationSectionDefinition)
  Section('REFERENCE SYSTEM', cls = ReferenceSystemSectionDefinition)
  Section('MAGNETISATION DIRECTION', [
    V('KMROT', int, 0),
    V('QMVEC', Array([0.,0.,0.])),
    V('DATA', Table({'MTET_Q' : float, 'MPHI_Q' : float }, numbering='IQ', free_header = True)),
    V('IT_DATA', Table({'IT': int, 'MTET_Q' : float, 'MPHI_T' : float, 'MGAM_T': float}, free_header = lambda x: '*' not in x),
      is_optional=True),
    ],
    cls = ArraySection('magnetisation_direction')
  )
  Section('MESH INFORMATION', cls = MeshInformationSectionDefinition)
  Section('OCCUPATION',   cls = OccupationSectionDefinition)
  Section('TYPES', cls = TypesSectionDefinition)

  return PotentialDefinition(sections)

potential_definition = fce()
""" Potential file format definition """

del fce
