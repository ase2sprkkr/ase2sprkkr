from ase.lattice.tetragonal import SimpleTetragonalFactory
from ase.visualize import view
from ase.spacegroup import get_spacegroup
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms
# Define a Perovskite Factory class
class PerovskiteFactory(SimpleTetragonalFactory):
    bravais_basis = [[0, 0, 0.0], [0.5, 0.5, 0.5],[0., 0.5, 0.5], [0.5, 0.5, 0],
                     [0.5, 0.0, 0.5]]
    element_basis = (0, 1, 2, 2, 2)

ABO3 = Perovskite = PerovskiteFactory()

# Generate the base STO cell
a0  = 3.905
STO = Perovskite(('Sr', 'Ti', 'O'),
                 latticeconstant={'a': a0, 'c/a': 1.},
                 size=(1,1,1))

#First we create new child of the atoms object which includes occupations
atoms=SPRKKRAtoms.promote_ase_atoms(STO)
atoms.sites[1].occupation.set({'Ti':0.99, 'Ni' : 0.01 })

calculator = SPRKKR(atoms=atoms,mpi=['mpirun','-np','4'])
calculator.input_parameters.set(NKTAB=50)
calculator.input_parameters.set(NL=3)
calculator.input_parameters.set(NE=32)
calculator.input_parameters.set(NITER=20)
calculator.input_parameters.SCF.MIX=0.01
calculator.input_parameters.ENERGY.ImE=0.0
calculator.input_parameters.ENERGY.GRID=5
out=calculator.calculate()
print(out.energy)
print(len(out.iterations))
print(out.iterations[-1]['error'])
print(out.last_iteration['moment'])

