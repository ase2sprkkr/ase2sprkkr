""" Here, the format of the potential_file is defined """

from .sections import \
    ScfInfoSectionDefinition, \
    GlobalSystemParameterDefinition, \
    LatticeSectionDefinition, \
    SitesSectionDefinition, \
    OccupationSectionDefinition, \
    ReferenceSystemSectionDefinition, \
    MeshInformationSectionDefinition, \
    TypesSectionDefinition, \
    PotentialSectionDefinition, \
    ChargeSectionDefinition, \
    HostMadelungPotentialSectionDefinition, \
    ChargeMomentsSectionDefinition, \
    MomentsSectionDefinition


def fce():
  from ...common.grammar_types import \
        Table, Array, Sequence, \
        Date, line_string
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
    V('SYSTEM', line_string, lambda x: 'System: {}'.format(x._get_root_container().atoms.symbols if x else '<UNKNOWN>') ),
    V('PACKAGE', line_string, 'SPR-KKR'),
    V('FORMAT', Sequence(int, Date(prefix='(', postfix=')'), names = ['VERSION', 'DATE']),
          default_value = [7, datetime.datetime(2007,5,21)]),
  ], name_in_grammar = False)

  Section('GLOBAL SYSTEM PARAMETER', cls = GlobalSystemParameterDefinition)
  Section('SCF-INFO', cls = ScfInfoSectionDefinition)
  Section('LATTICE', cls = LatticeSectionDefinition)
  Section('SITES', cls = SitesSectionDefinition)
  Section('OCCUPATION', cls = OccupationSectionDefinition)
  Section('REFERENCE SYSTEM', cls = ReferenceSystemSectionDefinition)
  Section('HOST MADELUNG POTENTIAL', cls = HostMadelungPotentialSectionDefinition)
  Section('CHARGE MOMENTS', cls = ChargeMomentsSectionDefinition)
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
  Section('OCCUPATION', cls = OccupationSectionDefinition)
  Section('TYPES', cls = TypesSectionDefinition)
  Section('POTENTIAL', cls = PotentialSectionDefinition)
  Section('CHARGE', cls = ChargeSectionDefinition)
  Section('MOMENTS', cls = MomentsSectionDefinition)

  return PotentialDefinition(sections)


potential_definition = fce()
""" Potential file format definition """

del fce
