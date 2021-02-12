from .sprkkr_io_types  import *
from .sprkkr_io_entry_def import \
    EntryDef as ED, \
    EntrySectionDef as ES, \
    EntryValueDef as V

def definition():

  CONTROL = ES[ 
      V('DATASET', str, required = True, help="Meaning of the parameter"),
      V('POTFIL', str, required=True, help="Potential file (see SPRKKR documentation for its format)"),
      V('KRWS', int)
  ]

  TAU = ES[
      V('BZINT', Keyword('POINTS', 'ANOTHER_OPTION'), required=true),
      V('NKTAB', int)
  ]

  ENERGY = ES[
      V('GRID', SetOf(int, length=1), required=True),
      V('NE', SetOf(int, length=1), required=True),
      V('ImE', float),
      V('EMIN', float),
  ]

  defs = {}
  defs['SCF'] = ED(
      CONTROL.copy(
        V('ADSI', Keyword('SCF'), required = True, help="Type of the computation -- do DFT selfconsistent cycle"),
      ),
      TAU,
      ENERGY.copy(
        ENERGY['GRID'].copy(default_value=1),
        ENERGY['NE'].copy(default_value=32)
      ),
      help = "SCF calculation"
  )

  defs['DOS'] = ED(
      CONTROL.copy(
        V('ADSI', Keyword('DOS'), required = True, help="Type of the computation"),
      ),
      TAU,
      ENERGY.copy(
          ENERGY['GRID'].copy(default_value=3),
      )
  )

  


