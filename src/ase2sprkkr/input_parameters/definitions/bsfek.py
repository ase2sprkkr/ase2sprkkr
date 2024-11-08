""" DOS task input parameters definition"""
from ...common.grammar_types import SetOf, Integer
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V
from ...common.generated_configuration_definitions import Length
from ...common.configuration_definitions import if_not_defined

input_parameters = lambda: InputParameters(
    'bsfek', [
          CONTROL('BSF'),
          TAU,
          TASK('BSF').copy([
            V('NK', 300, info="total number of k-points"),
            V('KPATH', Integer(min=1, max=5), info="Predefined path in k-space", description="""
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
""", is_optional=True),
            *if_not_defined('KPATH', [
                V('NKDIR', Length('KA', 'KE', default_values=[[0.,0.,0.],[1.,1.,1.,]]), info="Number of directions treated in k-spaces", is_required='Please, specify either TASK.KPATH or TASK.KA and TASK.KE'),
                V('KA', SetOf(float, length=3), is_repeated='NUMBERED', info="First k-vector segment in k-space in multiples of 2π/a and rectangular coordinates with * = 1, ...,NKDIR",
                is_required="Please, specify either TASK.KPATH or TASK.KA and TASK.KE"),
                V('KE', SetOf(float, length=3), is_repeated='NUMBERED', info="First k-vector segment in k-space in multiples of 2π/a and rectangular coordinates with * = 1, ...,NKDIR",
                is_required="Please, specify either TASK.KPATH or TASK.KA and TASK.KE"),
            ])
          ]),
          ENERGY(
            emin = (-0.2, 'The lowest E-value', None),
            emax = (-1., 'The highest E-value', None),
            defaults={
              'GRID': 3,
              'NE'  : 200,
              'ImE' : 0.001,
            }),
          CONTROL('BLOCHSF'),
          TAU,
          MODE,
          STRCONST,
          SITES
      ],
    executable='kkrgen',
    mpi=True,
    info="BSFEK - Bloch spectral functions in the E-K plane"
)
""" JXC -JXC task input parameters definition"""

# TODO - AKI scripts to generate KA/KE
