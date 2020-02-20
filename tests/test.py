# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:
import numpy as np
from sprkkr.calculator import SPRKKR

from ase.units import Bohr
from ase.build import bulk
from ase.lattice.compounds import L1_0
from ase.units import Bohr
from ase.lattice.tetragonal import SimpleTetragonalFactory

from sprkkr.misc import set_occupancy
from sprkkr.calcio import PotFile, InputFile


# Create a Perovskite, STO for example....
class PerovskiteFactory(SimpleTetragonalFactory):
    bravais_basis = [[0, 0, 0.5], [0.5, 0.5, 0],[0.5, 0, 0], [0, 0.5, 0], 
                     [0.5, 0.5, 0.5]]
    element_basis = (0, 1, 2, 2, 2)

ABO3 = Perovskite = PerovskiteFactory()

# Here is the base STO cell
a0 = 3.9
covera = 1
sto = Perovskite(('Sr', 'Ti', 'O'), latticeconstant={'a': a0, 'c/a': covera}, 
                 size=(1,1,1))

# ... and change the occupancy for the Ti site with some Ni ;-)
set_occupancy(sto, 1, 'Ti', 0.94)
set_occupancy(sto, 1, 'Ni', 0.06)


# create a SPRKKR calculator
calc = SPRKKR(label="calc/kkr")

# attach the atoms to the calculator object
calc.set_atoms(sto)

# Here is how to modify input file
calc.input.control_section.set(DATASET="SLA", ADSI="SLF", POTFIL="potfit_swap.pot_new")
calc.input.tau_section.set(nktab=50)
# Can also add a custom variable to a section, in that case, the value must be
# entered as a string
calc.input.tau_section.set(foo="2")
# can also add a custom section with custom variables...
calc.input.set('foobar', foo="3")

# launch kkrscf
calc.scf()
# This will write the input *.inp and *.pot file
# and run the kkrscf command


