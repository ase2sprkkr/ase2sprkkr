#!/usr/bin/env python
"""
This is a sctipt to visualise in_struct.inp files. Run it to see the doc.
"""
import numpy as np
import argparse
from pathlib import Path

if not __package__:
  __package__ = 'ase2sprkkr.tools'
import sys
sys.path.append(str(Path(__file__).resolve().parents[2]))


from ..potentials.potentials import Potential
from ..sprkkr.structure import structure_file_to_atoms
from ..ase.visualise import view

__author__='JM'

def visualise():
  #READ IN COMAND LINE PARAMETERS
  description='This is a sctipt to visualise in_struct.inp files. '
  description+='     You can either use ase visualisation tools or export it to cif file'
  parser = argparse.ArgumentParser(description=description)
  parser.add_argument('-i','--input', help='in_struct.inp file name (if not specified only pot file will be ploted)',required=False)
  parser.add_argument('-p','--pot', help='Filemane of scf Potential (always needed to be specified)',required=True)

  parser.add_argument('-o','--out',type=str,help='Output file name for structure (default structure.cif)', default='structure.cif',required=False)
  parser.add_argument('-f','--format',type=str,help='Output file fomrmat for structure (see ase allowed formats, default cif)', default='cif',required=False)
  parser.add_argument('-a','--ase',help='Use ase visualisation', action='store_true',required=False)
  parser.add_argument('-s','--scale-radii',help='Ase visualisation atomic radius', type=float,required=False, default=0.5)
  parser.add_argument('-v','--vac',dest='vacuum_height', type=float,help='Size of added vacuum in AA (default=10.0)',default=10.0,required=False)
  parser.add_argument('-b','--nbulk',type=int,help='Repetition of bulk unit (default=2)',default=2,required=False)

  args = parser.parse_args()
  ciffile=args.out
  cifpotfile=ciffile+'_pot'
  outformat=args.format
  ase_vis=args.ase
  structure_filename=args.input

  #read the potential
  potential = Potential.from_file(args.pot)
  pot_atoms = potential.atoms

  def ase_view(atoms):
      view(atoms, scale_radii=args.scale_radii)

  if (ase_vis):
      ase_view(pot_atoms)
  pot_atoms.write(cifpotfile, format = outformat)

  #From here on visualisation of in_structure.inp
  if not structure_filename:
     exit()
  # Extract data from structure file
  structure = structure_file_to_atoms(args.input, potential, n_bulk=args.nbulk, vacuum_height = args.vacuum_height)
  # Visualise the structure or write out to the file
  structure.write(ciffile, format = outformat)

  if (ase_vis):
      ase_view(structure)

if __name__ == "__main__":
  visualise()
