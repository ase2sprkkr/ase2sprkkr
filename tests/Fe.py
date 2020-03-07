# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:
import numpy as np
from ase.build import bulk
from sprkkr.calculator import SPRKKR

from ase.units import Bohr

from sprkkr.misc import set_occupancy
from sprkkr.calcio import PotFile, InputFile


# Create a Fe
Fe=bulk('Fe')

# create a SPRKKR calculator
calc = SPRKKR(label="Fe/Fe")

# attach the atoms to the calculator object
calc.set_atoms(Fe)

# Here is how to modify input file
#calc.input.control_section.set(DATASET="Fe", ADSI="SCF", POTFIL="Fe.pot")
#calc.input.tau_section.set(nktab=250)
# launch kkrscf
calc.get_potential_energy()
