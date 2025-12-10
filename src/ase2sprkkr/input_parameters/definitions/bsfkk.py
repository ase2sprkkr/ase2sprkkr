"""  bsfkk task input parameters definition"""
from ...common.grammar_types import SetOf
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V

input_parameters = lambda: InputParameters(
    'bsfkk', [
          CONTROL('BSF'),
          TAU,
          TASK('BSF', add=[
            V('NK1', int, info="number of k-points along k1", is_optional=True),
            V('NK2', int, info="number of k-points along k2", is_optional=True),
            V('K1', SetOf(float, length=3), is_optional=True, info="ﬁrst k-vector to span a two-dimensional region in k-space."),
            V('K2', SetOf(float, length=3), is_optional=True, info="second k-vector to span a two-dimensional region in k-space"),
          ]),
          ENERGY(
              emin = (0.7, 'the energy to compute the BSF', None),
              emax = 'emin',
              defaults={
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
    mpi=True,
    info="BSFKK - Bloch spectral functions in the K-K plane"
)
""" BSFKK - Bloch spectral functions in the K-K plane task input parameters definition"""
