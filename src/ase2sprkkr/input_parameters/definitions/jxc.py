"""JXC task input parameters definition"""
from .sections import TASK, CONTROL, TAU, ENERGY, SITES, STRCONST, MODE
from ..input_parameters_definitions import \
    InputParametersDefinition as InputParameters, \
    InputValueDefinition as V


def input_parameters():
    """ JXC -JXC task input parameters definition"""
    input_parameters = InputParameters(
        'jxc', [
              CONTROL('JXC'),
              TAU,
              MODE,
              STRCONST,
              ENERGY(),
              TASK('JXC').copy([
                  V('CLURAD', 2.2, info=
        """The radius of a sphere restricting the cluster
        of atoms for which the exchange coupling
        constants will be calculated"""
                    )
              ]),
              SITES
          ],
        executable='kkrgen',
        mpi=True,
        info="JXC"
    )
    input_parameters['MODE'].copy_member('MODE').warning_condition = lambda x: \
      "JXC task does not support SREL (scalar relativity without spin) or " \
      "NREL (no relativity at all) MODE. Please change SCF.MODE, or " \
      "the computation will fail." if x in ('SREL', 'NREL') else None
    input_parameters['CONTROL'].copy_member('NONMAG').warning_condition = lambda x: \
      "JXC task does not support non-magnetic computation. Please disable " \
      "CONTROL.NONMAG, or the computation will fail." if x else None

    return input_parameters
