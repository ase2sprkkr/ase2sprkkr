""" Here, the format of the potential_file is defined """

import datetime

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
    MagnetisationDirectionSectionDefinition, \
    MomentsSectionDefinition \


def _potential_header_default_value(_):
  return f'SPR-KKR potential file, created at {datetime.datetime.now()}'


def _potential_system_default_value(option):
  return 'System: {}'.format(option._get_root_container().atoms.symbols if option else '<UNKNOWN>')


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
    V('HEADER', line_string, _potential_header_default_value),
    Separator(),
    V('TITLE', line_string, 'Created by ASE-SPR-KKR wrapper'),
    V('SYSTEM', line_string, _potential_system_default_value),
    V('PACKAGE', line_string, 'SPR-KKR'),
    V('FORMAT', Sequence(int, Date(prefix='(', postfix=')'), names = ['VERSION', 'DATE']),
          default_value = [9, datetime.datetime(2019,1,18)]),
  ], name_in_grammar = False)

  Section('GLOBAL SYSTEM PARAMETER', cls = GlobalSystemParameterDefinition)
  Section('SCF-INFO', cls = ScfInfoSectionDefinition)
  Section('LATTICE', cls = LatticeSectionDefinition)
  Section('SITES', cls = SitesSectionDefinition)
  Section('OCCUPATION', cls = OccupationSectionDefinition)
  Section('REFERENCE SYSTEM', cls = ReferenceSystemSectionDefinition)
  Section('HOST MADELUNG POTENTIAL', cls = HostMadelungPotentialSectionDefinition)
  Section('CHARGE MOMENTS', cls = ChargeMomentsSectionDefinition)
  Section('MAGNETISATION DIRECTION', cls = MagnetisationDirectionSectionDefinition)
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
