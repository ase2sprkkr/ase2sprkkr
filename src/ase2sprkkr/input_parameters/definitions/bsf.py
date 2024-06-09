""" BSF task input parameters definition"""
from ...common.grammar_types  import *
from .sections import *
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V \

input_parameters = InputParameters(
      'bsf', [
      CONTROL('BLOCHSF'),
      TAU,
      MODE,
      ENERGY.copy([
        V('EMINEV', -5., info='Minimum of the energy window in eV with respect to the Fermi level'),
        V('EMAXEV', 0.5, info='Maximum of the energy window in eV with respect to the Fermi level'),
        ],
        remove = ['EMIN'],
        defaults = { 'GRID' : 3, 'NE' : 200 }
      ),
      SITES.copy(defaults = {'NL' : 4 }),
      TASK('BSF').copy([  V( 'NK', int,required=False, info = 'Number of K points along the path'), 
         V( 'KPATH', int, required=False, info='defines the k path for the band structure calculation'),   
         V( 'NKDIR', int, required=False, info='defines the number of segment used in your K path'),
         V( 'KA1', SetOf(float, length=3),required=False, info=''),
         V( 'KE1',SetOf(float, length=3),required=False,  info=''),
         V( 'KA2', SetOf(float, length=3),required=False, info=''),
         V( 'KE2', SetOf(float, length=3),required=False, info=''),
         V( 'KA3', SetOf(float, length=3),required=False, info=''),
         V( 'KE3', SetOf(float, length=3),required=False, info=''),
         V( 'KA4', SetOf(float, length=3),required=False, info=''),
         V( 'KE4', SetOf(float, length=3),required=False, info=''),
         V('KA',  Range(float), info='Scatering in momentum space ',required=False),
      V('K1', SetOf(float, length=3), info='Translating vector of the scatering in momentum space ',required=False),
      V('NK1', int, info='Number of momentum steps for the integration',required=False),
      V('K2',  SetOf(float, length=3), info='Translating vector 2 of the scatering in momentum space ',required=False),
      V('NK2', int, info='Number of momentum steps 2 for the integration',required=False),
      V('K3',  SetOf(float, length=3), info='Translating vector 3 of the scatering in momentum space ',required=False),
      V('NK3', int, info='Number of momentum steps 3 for the integration',required=False),
      V('K4',  SetOf(float, length=3), info='Translating vector 4 of the scatering in momentum space ',required=False),
      V('NK4', int, info='Number of momentum steps 4 for the integration',required=False),
         ]),
  ],
  executable='kkrgen',
  mpi=True
)
""" BSF task input parameters definition"""

from ...common.doc import process_input_parameters_definition
process_input_parameters_definition(__name__)
