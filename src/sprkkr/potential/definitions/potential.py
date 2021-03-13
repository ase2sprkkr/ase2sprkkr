from ..common.misc import OrderedDict
from .potential_defintions import \
      PotSectionDefinition, \
      PotValueDefinition as V,\
      Separator
from  ...common.grammar_types import *
from  ..lattice_section import LatticeSectionDefinition
from  ..sites_section import SitesSectionDefinition

sections = []
def Section(*args, cls=PotSectionDefinition, **kwargs):
   x=cls(name, values)
   sections.append(x)
   return x


Section('HEADER', [
  Separator(),
  V('HEADER', str),
  Separator(),
  V('TITLE', str),
  V('SYSTEM', str),
  V('PACKAGE', Keyword('SPRKKR')),
  V('FORMAT', Sequence(integer, Date(prefix='(', postfix=')'))),
], name_in_grammar = False)

Section('GLOBAL SYSTEM PARAMETER', [
  V('NQ', 1),
  V('NT', 2),
  V('NM', 1),
  V('IREL', 3)
])

Section('SCF-INFO', [
  V('INFO', str, None),
  V('SCFSTATUS', Keyword('START')),
  V('FULLPOT', False),
  V('BREITINT', False),
  V('NOMAG', False),
  V('ORBPOL', boolean, None),
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

Section('LATTICE',   cls = LatticeSectionDefinition)
Section('SITES',   cls = LatticeSectionDefinition)
