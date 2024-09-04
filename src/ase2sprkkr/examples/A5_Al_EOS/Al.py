""" Total energy SCF calculations as a function of Volume and corresponding
equation of state fiting. """

from ase.build import bulk
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase.io.trajectory import Trajectory
import numpy as np
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState


def main():
    al = bulk('Al', 'fcc', a=4.0)
    cell = al.get_cell()
    traj = Trajectory('Al.traj', 'w')

    # run calculation
    calc = SPRKKR(atoms=al,mpi=True)
    calc.input_parameters.set(NL=3)
    calc.input_parameters.set(NE=32)
    calc.input_parameters.SCF.MIX=0.20
    calc.input_parameters.ENERGY.ImE=0.0
    calc.input_parameters.ENERGY.GRID=5
    calc.input_parameters.SCF.FULLPOT = True
    calc.input_parameters.MODE.MODE='SREL'

    al.set_calculator(calc)

    for x in np.linspace(0.80, 1.20, 10):
       al.set_cell(cell * x, scale_atoms=True)
       al.get_potential_energy()
       traj.write(al)

    # Now lets plot the result and fit EOS
    configs = read('Al.traj@0:')  # read 5 configurations
    # Extract volumes and energies:
    volumes = [al.get_volume() for al in configs]
    energies = [al.get_potential_energy() for al in configs]
    eos = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    print(B / kJ * 1.0e24, 'GPa')
    eos.plot('Al-eos.png')


# Just run the script only when directly called from command line

if __name__ == "__main__":
    main()
