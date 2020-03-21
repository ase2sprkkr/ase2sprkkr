# coding: utf-8
import os
from msspec.calculator import MSSPEC
from msspec.utils import *

from ase.build  import bulk
from ase.visualize import view

from sprkkr.calculator import SPRKKR


a0 = 3.6 # The lattice parameter in angstroms

# Create the copper cubic cell
cluster = bulk('Cu', a=a0)
#
# FIRST  PERFORM SCF CALCULATIONS USING SPRKKR
#
calc_scf=SPRKKR(label="Cuscf/Cu",task='scf')
calc_scf.set_atoms(cluster)
calc_scf.get_potential_energy()

if calc_scf.converged:
   print("SCF Calculations converged after",calc_scf.niter,"iterations")
   conv_potfile=os.path.join(calc_scf.potfile+'_new')
else:
   raise "NOT CONVERGED"
#
# EXPORT POTENTIAL FOR PHAGEN
#
#change task and command
calc_scf.set_command('PHAGEN')
# Change output file
calc_scf.set_outfile('Cu_phagen.out')
#to change task we need to replace input file tempate
calc_scf.set_inpfile("Cu_phagen.inp")
calc_scf.input.control_section.set(DATASET="PHAGEN", ADSI="PHAGEN")
#set potetential file to converged potential
calc_scf.set_potfile(conv_potfile)
#run given task
calc_scf.phagen()

#
# Repeat the cell many times along x, y and z
cluster = cluster.repeat((8, 8, 8))
# Put the center of the structure at the origin
center_cluster(cluster)
# Keep atoms inside a given radius
cluster = cut_sphere(cluster, radius=a0 + .01)
# Keep only atoms below the plane z <= 0
cluster = cut_plane(cluster, z=0.1)

# Set the absorber (the deepest atom centered in the xy-plane)
cluster.absorber = get_atom_index(cluster, 0, 0, -a0)
# Create a calculator for the PhotoElectron Diffration
calc = MSSPEC(spectroscopy='PED')
# Set the cluster to use for the calculation
calc.set_atoms(cluster)

# Run the calculation
data = calc.get_theta_scan(level='2p3/2')

# Show the results
data.view()
#calc.shutdown()
