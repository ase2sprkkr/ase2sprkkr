#!/usr/bin/env python
import glob
import logging
import os
import sys
from msspec.calculator import MSSPEC
from msspec.utils import get_atom_index
from msspec.utils import hemispherical_cluster
from msspec.utils import SPRKKRPotential
from ase2sprkkr.sprkkr.calculator import SPRKKR

from ase.build import bulk



logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Create a copper cell
Cu = bulk('Cu')

# ########## SPRKKR part
if 'sprkkr' in sys.argv:
    # create a SPRKKR calculator
    calc = SPRKKR(atoms=Cu,mpi=['mpirun','-np','16'])

    # Here is how to modify input file
    

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
    #change task and command
    calc.input_parameters='PHAGEN'
#    calc.input_parameters.set(NL=3)
    calc.input_parameters.SCF.MIX=0.20
    calc.input_parameters.ENERGY.ImE=0.0
    calc.input_parameters.ENERGY.GRID=[5,3]
    calc.input_parameters.set(NE=32)
    calc.calculate(potential=out_scf.potential_filename)

# ######### MsSpec part
if 'msspec' in sys.argv:
    pot = SPRKKRPotential(Cu, "Cu_scf.pot_new", *glob.glob("*PHAGEN.pot"))

    nplanes = 3
    cluster = hemispherical_cluster(Cu, planes=nplanes,
                                    emitter_plane=nplanes-1)
    cluster.absorber = get_atom_index(cluster, 0, 0, 0)

    calc = MSSPEC(folder="calc")
    calc.set_atoms(cluster)
    calc.tmatrix_parameters.potential = pot

    data = calc.get_theta_scan(level='2p3/2')
    data.view()

if len(sys.argv) <= 1:
    print("Please specify either 'sprkkr', 'msspec' keywords or both "
          "of them on the command line")
                                                  
