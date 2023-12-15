""" Definitions of sections to be used in definitions of input parameters
(input files for SPR-KKR tasks) """

from ...common.grammar_types import DefKeyword, flag, energy, Integer, Array, Keyword

from ..input_parameters_definitions import \
      InputSectionDefinition as Section, \
      InputValueDefinition as V
from ...sprkkr.sprkkr_grammar_types import Site, AtomicType


def CONTROL(ADSI):
  """ Create the definition of the CONTROL section of the task input file.

  Parameters
  ----------
  ADSI: string
    the default value for the ADSI parameter of the resulting section

  Return
  ------
  CONTROL: InputSectionDefinition
  """

  return Section('CONTROL',[
      V('DATASET', str, 'case', required = True, result_is_visible=True,
                                info="The custom field for the description of the problem - the output files will have called 'DATASET.<ext>'."),
      V('ADSI', DefKeyword(ADSI), required = True, info="Type of the computation."),
      V('POTFIL', str, required=True, info="The potential file (see SPRKKR documentation for its format). It isn't necessary to set it, it will be set by the calculator."),
      V('KRWS', Integer(min=0, max=1), 1, info='If it is 0, RWS is taken from the potential file and scaled. If 1, RWS is calculated by scaling the muffin-tin radii by a common scaling factor. (This setting is forced in the case of FULLPOT.)'),
      V('KRMT', {
        '0' : 'RMT is taken from the potential file',
        '1' : 'RMT = min( x*RWS )',
        '2' : 'RMT = min( d_ij / 2 )',
        '3' : 'RMT from atomic charge density (=> KRWS=1)',
        '4' : 'RMT from atomic Hartree potential (=> KRWS=1)',
        '5' : 'RMT from total atomic potential (=> KRWS=1)',
        '6' : 'take average of 3 and 4 (=> KRWS=1)',
      }, default_value=None, info='It controls how the muffin-tin radii are calculated.',required=False),
      V('PRINT', Integer(min=0, max=5), default_value=0, info="Verbosity of the output (0-5). Do not affect the results in any way, just the amount of the printed output."),
      V('NONMAG', False, info="Set this flag, if it is known that the system considered is non-magnetic. This leads to a higher symmetry and a faster calculation. ")
  ])


class TauSection(Section):
    """ TAU Section of input parameters. """

    def set_from_atoms(self, section, atoms, io_data):
        """ NKTAB[xD] values has to be set according to the type of the problem."""
        for i in (section.NKTAB, section.NKTAB2D, section.NKTAB3D):
            i.clear_result()
        if atoms is not None:
           pbc = atoms.pbc.sum()
           # setting the dangerous value to the non-used options
           # prevents them to use its default values and so to be
           # writed to the output file
           if pbc == 3:
               section.NKTAB2D.result = section.NKTAB2D._create_dangerous_value(None)
               section.NKTAB3D.result = section.NKTAB3D._create_dangerous_value(None)
           else:
               section.NKTAB.result = section.NKTAB._create_dangerous_value(None)


def _nktab_value(option):
    return option._container.NKTAB()


TAU = TauSection('TAU',[
      V('BZINT', DefKeyword({'POINTS' : 'special points method',
                             'WEYL' : 'Weyl method'},
                             description=
"""
The Weyl method (BZINT=WEYL) is a point sampling method using more or less ran-
dom points. The number of k-points used for the integration varies quadratically be-
tween 0.0 and ImE according to the imaginary part of the energy.

The special point method (BZINT=POINTS) uses a regular k-point grid with NKTAB
points. It is the standard method and gives a good compromise concerning accuracy
and efficiency. For BZINT=POINTS the parameter NKTAB will be adjusted to allow a
regular mesh.
"""), required=True,
                             info='The mode of BZ-integration used for calculation of the scattering '
                                  ' path operator τ'),
      V('NKTAB', 250, info='Number of points for the special points method', is_optional=True,
        description='For 2D problem, it serves only as default value for NKTABxD.'),
      V('NKTAB2D', int, _nktab_value, info='Number of points for the special points method for 2D region of 2D problem',
        is_optional=True, description='If it is not specified, NKTAB is used.'),
      V('NKTAB3D', int, _nktab_value, info='Number of points for the special points method for 3D region of 2D problem',
        is_optional=True, description='If it is not specified, NKTAB is used'),
      V('NKMIN', 300, info='Minimal number of k-points used for Weyl integration', write_condition= lambda o: o._container.BZINT() == 'WEYL'),
      V('NKMAX', 500, info='Maximal number of k-points used for Weyl integration', write_condition= lambda o: o._container.BZINT() == 'WEYL'),
      # expert
      V('CLUSTER', flag, expert=False, info="""Do cluster type calculation.""", description=
        "Cluster type calculation calculate τ by inverting the real space KKR matrix. "
        "Specify cluster center using IQCNTR or ITCNTR and its size using NSHLCLU or CLURAD."
      ),
      V('NSHLCLU', int, is_expert=True, is_optional=True, info="Number of atomic shells around the central atom of a cluster"),
      V('CLURAD', int, is_expert=True, is_optional=True, info="Radius of the cluster in multiples of ALAT."),
      V('IQCNTR', Site.I, is_expert=True, is_optional=True, info="The center of the cluster is set at the site position with number IQCNTR of the specified basis."),
      V('ITCNTR', AtomicType.I, is_expert=True, is_optional=True, info="The center of the cluster is set at one of the site positions that is occupied by the atomic type ITCNTR."),
      V('NLOUT', expert=3, info="The calculated τ -matrix is printed up to lmax=NLOUT."),
      V('MOL', expert=False, info="Cluster type calculation but for a molecular system. The system is specified as for CLUSTER.")
  ])
"""The definition of the TAU section of the task input file """

ENERGY = Section('ENERGY',[
      V('GRID', [5], required=True),
      V('NE', [32], required=True, info='Number of points in energy-mesh'),
      V('ImE', energy, 0.0),
      V('EMIN', -0.2, info='The real part of the lowest E-value'),
  ])
"""The definition of the ENERGY section of the task input file """

SCF = Section('SCF', [
      V('NITER', 200, info='Maximal number of iterations of the SCF cycle'),
      V('MIX', 0.2, info='Mixing parameter'),
      V('MIXOP', float, required=False),
      V('VXC', DefKeyword({
        'VWN' : 'Vosko, Wilk, Nusair',
        'MJW' : 'Janak, Williams, Moruzzigit g',
        'VBH' : 'von Barth, Hedin',
        'PBE' : 'Perdew, Burke, Ernzendorfer GGA',
        'PW92' : 'Perdew Wang',
        'EV-GGA' : 'Engel and Vosko GGA',
        'BJ' : 'Becke-Johnson',
        'MBJ' : 'modified Becke-Johnson'

        }), info='parametrisation of the exchange-correlation potential'),
      V('ALG', DefKeyword({
          'BROYDEN2': 'Broyden’s second method',
          'TCHEBY': 'Tchebychev'
        }), info='Mixing algorithm'),
      V('EFGUESS', float, required=False, info='Skip the Fermi energy search in the beginning.'),
      V('TOL', 0.00001, info='Tolerance threshold for the mixing algorithm'),
      V('ISTBRY', 1, info='Start Broyden after ISTBRY iterations'),
      V('FULLPOT', False, info='Non-spherical callculation (full-potential) instead of ASA'),
      V('ITDEPT', 40, info='Iteration depth for Broyden algorithm (length of used history)'),
      V('QION', Array(float), required=False, info='Guess for the ionic charges Qt for atomic types'),
      V('MSPIN', Array(float), required=False, info='Guess for the magnetic moment μ_{spin,t} for atomic types'),
      V('USEVMATT', False, info='Set up the starting potential using the original Mattheiss'
                                 'construction for the potential V instead of the charge density')
  ])
"""The definition of the SCF section of the task input file """

STRCONST = Section('STRCONST', [
        V('ETA', float, is_optional=True, info='Ewald parameter'),
        V('RMAX', float, is_optional=True, info='Convergency radius in real space'),
        V('GMAX', float, is_optional=True, info='Convergency radius in reciprocal space'),
      ], is_expert=True, is_optional=True, info=
      """The calculation of the ~k-dependent KKR structure constant matrix G(~k, E) is controlled by
three convergence parameters. ETA determines the relative weight of the real and reciprocal
space lattice sums, that are determined by the convergence radii RMAX and GMAX, respec-
tively. These convergence parameters have to be optimised anew if the lattice structure, the
lattice parameter or the energy or ~k-range used is changed. This is done by the program if
no values are applied via the input file. In some cases, in particular if one works at high
energies, it might be necessery to set the convergence by hand. For this purpose one can start
from the values set by kkrgen or kkrscf (see the output file)."""
)
"""The definition of the STRCONST section of the task input file """

SITES = Section('SITES', [
      V('NL', [3], info='Angula momentum cutoff (the first discarded l-space)', description=
""" The KKR-method is a minimum basis set method. This means that the angular-momentum
expansion can be chosen according to the atomic properties of the atomic types. For tran-
sition metals it is therefore normally sufficient to have a maximum l-value of 2 (NL = 3).
For systems with many atoms per unit cell it is in principle possible to set the l-expansion
according to the atom types on the lattice sites. For US having the NaCl-structure one could
choose NL = 4 for the U-site and NL = 2 for the S-site. At the moment this possibility, that
would save storage and computer time, is not supported by all subroutines. For this rea-
son a common l-expansion cutoff is used, that is fixed by the highest that occurs. For US
this implies that NL = 4 is used for all sites.""")
  ])
"""The definition of the SITES section of the task input file """


CPA = Section('CPA', [
    V('NITER', 20, info="Maximum number of CPA iterations"),
    V('TOL', 0.0001, info="Threshold for stopping CPA-cycle"),
  ], is_expert=True, is_optional=True,
  info="""For a system with substitutional disorder, the CPA is used. The listed variables control the CPA cycle """,
)
"""CPA section definition"""


MODE = Section('MODE', [
    V('MODE',
      Keyword({
        'NREL' : "work in the nonrelativistic mode",
        'SREL' : "work in the scalar-relativistic mode",
        'SP-SREL' : "work in the spin-polarized scalar-relativistic mode"}, aliases= {'FREL' : None }
      ), None,
                required=False, name_in_grammar=False,
                info='Using this option you can switch on the spin polarization and relativistic mode. If it''s not set (or set to ''FREL''), the ''full'' relativity mode is used.'),
    V('LLOYD', False, info='Use LLoyd formula for scattering operator. It can improve the accuracy of the Fermi energy.'),
    V('MDIR', Array(float, length=3), [1.,0.,0.], required=False, info="Common magnetisation direction vector with x, y and z in Cartesian coordinates. The normalisation is arbitrary.", is_numbered_array=True, is_always_added=False ),
    V('C', 1.0, info='Scale the speed of light for a given atom type.', is_numbered_array=True, required=False, is_always_added=False),
    V('SOC', 1.0, info='Scale the strength of the spin-orbit coupling for atom type.', is_numbered_array=True, required=False, is_always_added=False),
  ], is_expert=True, is_optional=True, info=
      """This section contains options that describe, how to consider relativity and/or spin. If the MODE is not specified otherwise the programs of the SPRKKR-package assume that a magnetic system should be treated in a fully relativistic way. By setting the parameter SP-SREL a slightly faster scalar relativistic calculation can be done instead for a magnetic system.""",
)
"""MODE Section definition"""

SPEC=Section('SPEC',[
    V('FEGFINAL', flag, required=False)
  ], is_expert=True, is_optional=True, info='Sets the final state as a free electron gas final state. By default it is TRLEED final state which is time reverse low energy electron diffraction final state.'
)


def TASK(task, add = []):
  """ Create the definition of the CONTROL section of the task input file.

  Parameters
  ----------
  TASK: string
    the default value for the TASK parameter of the resulting section

  Return
  ------
  TASK: InputSectionDefinition
  """
  if isinstance(task, str):
     task = DefKeyword(task)

  return Section('TASK', [
    V('TASK', task, name_in_grammar=False)
  ] + add)


__all__ = [
    'CONTROL', 'TAU', 'ENERGY', 'SCF', 'SITES', 'TASK', 'STRCONST', 'CPA', 'MODE'
]
