""" DOS task input parameters definition"""
from ...common.grammar_types import SetOf
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V
from ...common.doc import process_input_parameters_definition

input_parameters = InputParameters(
    'bsfkk', [
          CONTROL('BSF'),
          TAU,
          TASK('BSF').copy([
            V('NK', 300, info="total number of k-points"),
            V('NK1', int, info="number of k-points along k1", is_optional=True),
            V('NK2', int, info="number of k-points along k2", is_optional=True),
            V('KA', SetOf(float, length=3), default_value=[0.,0.,0.], info="Shift in the k-space."),
            V('K1', SetOf(float, length=3), is_optional=True, info="ﬁrst k-vector to span a two-dimensional region in k-space."),
            V('K2', SetOf(float, length=3), is_optional=True, info="second k-vector to span a two-dimensional region in k-space"),
          ]),
          ENERGY.copy([
              V('EMAX', 1., info="highest E-value"),
            ], defaults={
              'EMIN': 0.7,   # TODO - fermi energy from computation
              'ImE' : 0.001,
              'GRID': 3,
              'NE'  : 1
          }),
          CONTROL('BLOCHSF'),
          TAU,
          MODE,
          STRCONST,
          SITES
      ],
    executable='kkrgen',
    mpi=True
)
""" JXC -JXC task input parameters definition"""

process_input_parameters_definition(__name__)

# TODO - AKI scripts to generate KA/KE