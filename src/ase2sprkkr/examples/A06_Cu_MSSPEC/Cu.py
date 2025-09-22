#!/usr/bin/env python
"""
Simple SCF calculations for fcc Cu in combination with msspec code to generate
photoelectron diffraction.

Call with ``msspec`` and/or ``sprkkr`` commandline argument.
"""

import glob
import sys
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase.build import bulk


def main(args):
    # Create a copper cell
    Cu = bulk('Cu')

    # ########## SPRKKR part
    if 'sprkkr' in args:
        # create a SPRKKR calculator
        calc = SPRKKR(atoms=Cu,mpi=True)

        # launch kkrscf
        calc.input_parameters.set(NL=3)
        calc.input_parameters.SCF.MIX=0.20
        calc.input_parameters.ENERGY.ImE=0.0
        calc.input_parameters.ENERGY.GRID=[5,3]
        calc.input_parameters.set(NE=32)
        out_scf=calc.calculate()
        #
        # EXPORT POTENTIAL FOR PHAGEN
        #

        # Now use the calculator with the newly converged potential
        calc = out_scf.calculator
        calc.input_parameters='PHAGEN'
        calc.input_parameters.SCF.MIX=0.20
        calc.input_parameters.ENERGY.ImE=0.0
        calc.input_parameters.ENERGY.GRID=[5,3]
        calc.input_parameters.set(NE=32)
        calc.calculate(potential=out_scf.potential_filename)

    # ######### MsSpec part
    if 'msspec' in args:

        from msspec.calculator import MSSPEC
        from msspec.utils import get_atom_index, hemispherical_cluster, \
                                 SPRKKRPotential
        pot = SPRKKRPotential(Cu, "Cu_scf.pot_new", *glob.glob("*PHAGEN.pot"))

        nplanes = 3
        cluster = hemispherical_cluster(Cu, planes=nplanes,
                                        emitter_plane=nplanes - 1)
        cluster.absorber = get_atom_index(cluster, 0, 0, 0)

        calc = MSSPEC(folder="calc")
        calc.set_atoms(cluster)
        calc.tmatrix_parameters.potential = pot

        data = calc.get_theta_scan(level='2p3/2')
        data.view()

    if len(args) <= 1:
        print("Please specify either 'sprkkr', 'msspec' keywords or both "
              "of them on the command line")


# Just run the script only when directly called from command line

if __name__ == "__main__":
    main(sys.argv)
