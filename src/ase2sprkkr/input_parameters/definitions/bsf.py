""" DOS task input parameters definition"""
from ...common.grammar_types import Array, Integer
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE, SCF
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V
from ...common.doc import process_input_parameters_definition

input_parameters = InputParameters(
    'bsf', [
          CONTROL('BSF'),
          TAU,
          TASK('BLOCHSF'),
          ENERGY.copy([
            V('EMAX', float, info="highest E-value"),
            V('NK', 51, info="total number of k-points"),
            V('KPATH', Integer(min=1, max=5), 1, info="Predefined path in k-space", description="""
Bravais-lattice KPATH path
==========================
orb 1  Γ-Σ-X-G-U-A-Z-Λ-Γ-∆-Y-H-T-B-Z
       + X-D-S-C-Y + U-P-R-E-T + S-Q-T
    2  Γ-Σ-X-G-U-A-Z-Λ-Γ-∆-Y-H-T-B-Z
    3  Γ-Σ-X-G-U-A-Z-Λ-Γ
    4  Γ-∆-Y-H-T-B-Z
hex 1  Γ-Σ-M-T’-K-T-Γ-∆-A-R-L-S’-H-S-A
       + M-U-L + K-P-H
    2  Γ-Σ-M-T’-K-T-Γ-∆-A-R-L-S’-H-S-A
    3  Γ-Σ-M-T’-K-T-Γ-∆-A
    4  Γ-Σ-M
    5  K-T-Γ
sc  1  Γ-∆-X-Y-M-V-R-Λ-Γ-Σ-M
    2  Γ-∆-X-Y-M-V-R-Λ-Γ
    3  Γ-∆-X-Y-M-V-R
    4  Γ-∆-X-Y-M
fcc 1  X-∆-Γ-Λ-L-Q-W-N-K-Σ-Γ
       + L-M-U-S-X-Z-W-D-U
    2  X-∆-Γ-Λ-L-Q-W-N-K-Σ-Γ
    3  X-∆-Γ-Λ-L
    4  Γ-∆-X
    5  Γ-Λ-L
bcc 1  Γ-D-H-G-N-Σ-Γ-Λ-P-F-H + N-D-P
    2  Γ-D-H-G-N-Σ-Γ-Λ-P-F-H
    3  Γ-D-H-G-N-Σ-Γ-Λ-P
    4  Γ-D-H-G-N-Σ-Γ
    5  Γ-D-H
"""),
            V('NKDIR', 1, info="Number of directions treated in k-spaces"),
            V('KA', Array(float, length=3), [0.,0.,0.], info="First k-vector segment in k-space in multiples of 2π/a and rectangular coordinates with * = 1, ...,NKDIR"),
            V('KE', Array(float, length=3), [1.,0.,0.], info="First k-vector segment in k-space in multiples of 2π/a and rectangular coordinates with * = 1, ...,NKDIR"),
            V('NK1', int, info="number of k-points along k1", is_optional=True),
            V('NK2', int, info="number of k-points along k2", is_optional=True),
            V('K1', Array(float, length=3), [1.,0.,0.], info="ﬁrst k-vector to span a two-dimensional region in k-space."),
            V('K2', Array(float, length=3), [1.,0.,0.], info="second k-vector to span a two-dimensional region in k-space"),
          ], defaults={
            'ImE' : 0.1
          }),

          CONTROL('BLOCHiSF'),
          TAU,
          MODE,
          STRCONST,
          SCF,
          SITES
      ],
    executable='kkrgen',
    mpi=True
)
""" JXC -JXC task input parameters definition"""

process_input_parameters_definition(__name__)
