import os
from ase.build import bulk
from ase.calculators.calculator import kpts2mp
from sprkkr.calculator import SPRKKR
from ase.io.trajectory import Trajectory
import numpy as np


al = bulk('Al', 'fcc', a=4.0)
cell = al.get_cell()
traj = Trajectory('Al.traj', 'w')

# run calculation
calc = SPRKKR(label='Al',task='scf')
al.set_calculator(calc)

for x in np.linspace(0.80, 1.20, 10):
   al.set_cell(cell * x, scale_atoms=True)
   al.get_potential_energy()
   traj.write(al)


