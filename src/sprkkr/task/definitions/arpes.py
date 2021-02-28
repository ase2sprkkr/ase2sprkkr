from .option_types  import *
from .sections import *
from .configuration_definitions import \
    TaskDefinition as Task, \
    OptionDefinition as O

taks = Task(
    'arpes', [
    CONTROL('ARPES'),
    TAU,
    ENERGY.copy([
        O('EMINEV', -10.),
        O('EMAXEV', -10),
        O('EWORKEV', 4.2),
        O('IMV_INI_EV', 0.05),
        O('IMV_FIN_EV', 5.),
      ],
      remove = ['EMIN'],
      defaults = { 'GRID' : 1 }
    ),
    SITES.copy(defaults = {'NL' : 4 }),
    TASK('ARPES').copy([
      O('STRVER', 0),
      O('IQ_AT_SURF', 2),
      O('MILLER_HKL', SetOf(Integer, length=0), [0,0,1]),
      O('NTMP', 1),
      O('TMPMIN', 11.),
      O('CTMPMAX', 11.),
      O('CTMPMAX', 11.),
      O('VIBRA', Flag, True),
      O('CNVIBRA', 14),
    ]),
    Section('SPEC_PH', [
      O('THETA', 45.),
      O('PHI', 0.),
      O('POL_P', 'P'),
      O('EPHOT', 6675),
    ]),
    Section('SPEC_STR', [
      O('N_LAYDBL', SetOf(int), [10,10]),
      O('NLAT_G_VEC', 57),
      O('N_LAYER', 120),
      O('SURF_BAR', SetOf(float), [0.25,0.25]),
      O('TRANSP_BAR', Flag, False)
    ])
  ])

