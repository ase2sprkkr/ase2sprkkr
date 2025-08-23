# Script demonstrating how to set up a disordered local moment (DLM) calculation
# using bcc Fe as the test case.
#
# Such calculations are used to describe the paramagnetic state of materials,
# where local magnetic moments on atoms persist, but where the net magnetisation
# of the material is zero on account of thermally induced spin fluctuations. The 
# CPA is used to average over possible orientations of magnetic moment. Above the Curie
# temperature, all orientations of magnetic moment are equally probably, and the
# probability distribution of orientations of magnetic moments is uniform on the sphere. 
# In the absence of spin-orbit coupling, the CPA condition for a full average (using 
# an angular integral) over the sphere with uniform probability distribution can be shown 
# to be mathematically equivalent to the CPA condition for a 50:50 Ising 'alloy' of 'up' 
# and 'down' magnetic moments. (Note that this equality does not strictly hold in the 
# case where spin-orbit effects are present, as in a fully relativistic calculation.)
#
# A DLM calculation should be iterated to self-consistency to check whether the local
# moments 'collapse' (as frequently happens for elements like Cr, Ni) or whether they
# are self-sustaining (as is typical for Fe, Co).
#
# Once you have run this script and iterated to self-consistency, if you run
# `tail Fe.pot_new`
# you should find that, although the net magnetisation of the calculation is zero,
# the two Fe atoms have equal and opposite magnetic moments, each of magnitude
# approximately 2.16 \mu_B. This confirms that the local moments are self-sustaining,
# and that bcc Fe is a 'good' local moment system.
#
# Some appropriate references for the disordered local moment (DLM) picture are:
# - A. J. Pindor et al., J. Phys. F: Met. Phys. 13, 979 (1983)
# - J. B. Staunton et al., J. Magn. Magn. Mater. 45, 15-22 (1984)
# - B. L. Gyorffy et al., J. Phys. F: Met. Phys. 15, 1337 (1985)
#
# Author: Christopher D. Woodgate, University of Bristol, 2025
#
# Email: christopher.woodgate@bristol.ac.uk
def main():

    import numpy as np
    from ase import Atoms
    from ase.build import bulk
    from ase2sprkkr.sprkkr.calculator import SPRKKR
    from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms

    # BCC Fe has a_lat = 2.87 \AA
    a_lat = 2.87
    atom = bulk('Fe', 'bcc', a=a_lat)

    # Now set up the 'alloy'
    atoms=SPRKKRAtoms.promote_ase_atoms(atom)

    # We cannot use a Python dictionary here, as we need an alloy of two elements
    # with the same chemical symbol and this behaviour is not permitted in a Python dictionary
    # Instead, we use a list.
    atoms.sites[0].occupation.set([('Fe', 0.5), ('Fe',0.5)])
    
    # Now we set up the ASE calculator
    calculator = SPRKKR(atoms=atoms,mpi=['mpirun','-np','4'])

    # NOTE: we need to specify initial magnetic moments, up and down, on the two Fe atoms
    # NOTE: These are in units of \mu_B
    calculator.input_parameters.SCF.MSPIN = [2.3, -2.3]

    # And other necessary parameters
    calculator.input_parameters.set(NKTAB=500) # Control k-point sampling
    calculator.input_parameters.set(VXC='VWN') # VWN LDA XC functional
    calculator.input_parameters.set(MODE='SP-SREL') # Scalar relativistic mode for the DLM physics
    calculator.input_parameters.set(NL=4) # l_max=3 means 4 l-channels
    calculator.input_parameters.set(NE=32) # 32-point semicircular energy contour
    calculator.input_parameters.ENERGY.GRID = 5 # Semicircular energy contour
    calculator.input_parameters.SCF.MIX = 0.05 # Gentle mixing paramter to ensure convergence
    calculator.input_parameters.SCF.TOL = 1e-7 # Tight tolerance on SCF cycle
    calculator.input_parameters.SCF.NITER = 1000 # Allow for lots of iterations if needed
    calculator.input_parameters.ENERGY.ImE = 0.005 # Distance of closest approach to real axis
    calculator.input_parameters.CPA.NITER = 100 # Allow for a lot of CPA iterations if needed
    calculator.input_parameters.CPA.TOL = 1e-8 # Tight CPA tol to go with tight SCF tol
    calculator.input_parameters.SCF.USEVMATT = True # I like to use the Mattheiss construction for the potential

    # Run the calculator
    out=calculator.calculate()

    print(out.energy)
    print(out.last_iteration['energy']['EF'].to_dict())
    print(out.iterations[-1]['error']())
    print(out.last_iteration['moment'].to_dict())

# If we run this script, execute 'main'
if __name__ == "__main__":
    main()
