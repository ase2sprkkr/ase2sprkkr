""" XAS task input parameters definition"""
from ...common.grammar_types import SetOf, DefKeyword, Keyword, Range, Array, flag
from .sections import CONTROL,TASK, TAU, ENERGY, SITES, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V, \
    InputSectionDefinition as Section


def input_parameters():
    """ XAS -xas task input parameters definition"""
    input_parameters = InputParameters('xas', [
              CONTROL('XAS'),
              TAU,
              MODE,
              ENERGY(defaults={
                'EMAX':4.0,
                'ImE' : 0.01,
                'GRID': 6,
                'NE'  : 180}),
              TASK('XAS').copy([ V('IT', int,1, is_required=True, info= """atom type IT"""),
                  V('CL', str,'2P',is_required=True, info="""initial core level shell"""),
                  V('MECHECK', flag, False, is_optional=True ),
                  V('OUTPUT',DefKeyword({'MBARN':'output of absorption coefficient µas µ atom = µ Vuc in [Mbarn] ',
                    'SIGMA':'SIGMA : output as absorptive part of optical conductivity σ = µ c/4π in [10E15/s ]'}),
                    is_expert=True, is_optional=True, info='write extra output'),
                  V('FRAMETET', float, is_expert=True, is_optional=True,
                  info='Polar angle θ (FRAMETET) defining the orientation of the electric field vector of the incident light with respect to the material surface normal. Default is FRAMETET = 0, meaning the field lies along the surface normal.'),
                   V('FRAMEPHI', float, is_expert=True, is_optional=True, 
                   info='Azimuthal angle φ (FRAMEPHI) defining the in-plane rotation of the electric field vector of the incident light relative to the surface reference axis. Default is FRAMEPHI = 0, corresponding to alignment with the x-axis of the surface frame.'),
              ]),
              SITES
          ],
        executable='kkrgen',
        mpi=True,
        info="Calculates X-ray absorption spectra"
    )
    return input_parameters
