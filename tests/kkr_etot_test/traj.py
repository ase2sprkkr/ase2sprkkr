from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
configs = read('Al.traj@0:')  # read 5 configurations
# Extract volumes and energies:
volumes = [al.get_volume() for al in configs]
energies = [al.get_potential_energy() for al in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('Al-eos.png')
