# coding: utf-8
# vim: set et sw=4 ts=4 sts=4 nu ai cc=80 fdm=indent mouse=a:
import numpy as np
from sprkkr.calculator import SPRKKR

from ase.units import Bohr
from ase.build import bulk
from ase.lattice.compounds import L1_0
from ase.units import Bohr
from ase.lattice.tetragonal import SimpleTetragonalFactory

from sprkkr.misc import set_occupancy,get_occupancy
from sprkkr.calcio import PotFile, InputFile
from ase.spacegroup import Spacegroup,get_spacegroup
from ase.build  import cut


def primitive_from_conventional_cell(atoms, spacegroup=1, setting=1):
    """Returns primitive cell given an Atoms object for a conventional
    cell and it's spacegroup."""
    from ase.spacegroup import Spacegroup
    from ase.build  import cut
    sg = Spacegroup(spacegroup, setting)
    prim_cell = sg.scaled_primitive_cell  # Check if we need to transpose
    return cut(atoms, a=prim_cell[0], b=prim_cell[1], c=prim_cell[2])

# Create a Perovskite, STO for example....
class PerovskiteFactory(SimpleTetragonalFactory):
    bravais_basis = [[0.5, 0.5, 0.0], [0.5, 0.5, 0.5],[0.0, 0.0, 0], [0, 0.5, 0.5], 
                     [0.5, 0.0, 0.5]]
    element_basis = (2, 1, 0, 2, 2)

ABO3 = Perovskite = PerovskiteFactory()

# Here is the base STO cell
a0 = 3.9
covera = 1
sto = Perovskite(('Sr', 'Ti', 'O'), latticeconstant={'a': a0, 'c/a': covera}, 
                 size=(1,1,1))
# ... and change the occupancy for the Ti site with some Ni ;-)
set_occupancy(sto, 0, 'O',  0.94)
set_occupancy(sto, 0, 'Ni', 0.06)
#sto2=primitive_from_conventional_cell(sto)
#sg = get_spacegroup(sto)
#siteso = sto.get_scaled_positions()
#sites = sto.get_scaled_positions()
#sites=sg.unique_sites(sites)
#print(sites,"\n\n\n")
#a=[]
#for i in range(len(sites)):
# s,k=sg.equivalent_sites(sites[i])
# a.append(s) 
#tol=1e-5
#k=0
#typ=np.zeros(shape=(len(sites),len(siteso)),dtype=int)
#for site in a:
# for i in range(len(site)):
#  for j in range(len(siteso)):
#     dif=np.abs(site[i]-siteso[j]) <= [tol,tol,tol]
#     if dif.all():
#        typ[k][j]=1
# k=k+1
#
#occupancy = get_occupancy(sto)
#print(occupancy)
#print(a)
#for a in sto:
#    print(a)
#ntypes=typ.shape[0]
#print(typ,typ.shape,ntypes)
#exit()
#
#
#
# create a SPRKKR calculator
calc = SPRKKR(label="calc/kkr")

# attach the atoms to the calculator object
calc.set_atoms(sto)

# Here is how to modify input file
calc.input.control_section.set(DATASET="SLA", ADSI="SLF", POTFIL="kkr.pot")
calc.input.tau_section.set(nktab=50)
# Can also add a custom variable to a section, in that case, the value must be
# entered as a string
calc.input.tau_section.set(foo="2")
# can also add a custom section with custom variables...
calc.input.set('foobar', foo="3")

# launch kkrscf
#calc.scf()
calc.get_potential_energy()
# This will write the input *.inp and *.pot file
# and run the kkrscf command


