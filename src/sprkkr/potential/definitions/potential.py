def fce():
  from ...common.misc import OrderedDict
  from  ...common.grammar_types import DefKeyword, Table, Array, Sequence, integer, Date, boolean, line_string
  from ..potential_definitions import \
        PotSectionDefinition, \
        ASEArraySectionDefinition, \
        PotValueDefinition as V,\
        PotentialDefinition, \
        Separator
  from  .lattice import LatticeSectionDefinition
  from  .sites import SitesSectionDefinition
  from  .occupation import OccupationSectionDefinition
  from  .reference_system import ReferenceSystemSectionDefinition
  from  .types import TypesSectionDefinition
  from  .mesh_information import MeshInformationSectionDefinition
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
    V('HEADER', line_string, lambda x: f'SPRKKR potential file, created at {datetime.datetime.now()}'),
    Separator(),
    V('TITLE', 'SPR-KKR '),
    V('SYSTEM', line_string, lambda x: 'System: {}'.format(x._get_root_container().atoms.symbols) ),
    V('PACKAGE', DefKeyword('SPRKKR')),
    V('FORMAT', Sequence(int, Date(prefix='(', postfix=')')),
          default_value = [7, datetime.datetime(2007,5,21)]),
  ], name_in_grammar = False)

  Section('GLOBAL SYSTEM PARAMETER', [
    V('NQ', 1),
    V('NT', 2),
    V('NM', 1),
    V('IREL', 3)
  ])

  Section('SCF-INFO', [
    V('INFO', str, 'NONE'),
    V('SCFSTATUS', DefKeyword('START')),
    V('FULLPOT', False),
    V('BREITINT', False),
    V('NOMAG', False),
    V('ORBPOL', str, 'NONE'),
    V('EXTFIELD', False),
    V('BLCOUPL', False),
    V('BEXT', 0.0),
    V('SEMICORE', False),
    V('LLOYD', False),
    V('SCF-ITER', 0),
    V('SCF-MIX', 0.2),
    V('SCF-TOL', 1e-5),
    V('RMSAVV', 999999.),
    V('RMSAVB', 999999.),
    V('EF', 999999.),
    V('VMTZ', 0.7),
  ])

  Section('LATTICE', cls = LatticeSectionDefinition)
  Section('SITES',   cls = SitesSectionDefinition)
  Section('OCCUPATION',   cls = OccupationSectionDefinition)
  Section('REFERENCE SYSTEM', cls = ReferenceSystemSectionDefinition)
  Section('MAGNETISATION DIRECTION', [
    V('KMROT', int, 0),
    V('QMVEC', Array([0.,0.,0.])),
    V('DATA', Table({'QMTET' : float, 'QMPHI' : float }, numbering='IQ')),
    ],
    cls = ArraySection('magnetisation_direction')
  )
  Section('MESH INFORMATION', cls = MeshInformationSectionDefinition)
  Section('OCCUPATION',   cls = OccupationSectionDefinition)
  Section('TYPES', cls = TypesSectionDefinition)

  return PotentialDefinition(sections)

potential_definition = fce()
del fce
