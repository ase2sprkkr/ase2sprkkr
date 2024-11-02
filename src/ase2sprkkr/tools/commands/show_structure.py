#!/usr/bin/env python
"""
This is a sctipt to visualise potential or in_struct.inp files. Run it to see the doc.
"""
from pathlib import Path

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
import sys
sys.path.append(str(Path(__file__).resolve().parents[3]))

help='Visualise potential (and possibly in_struct.inp) files.'
description='You can either use ase visualisation tools or export it to cif file'


def parser(parser):
  parser.add_argument('potential_file', help='SPRKKR potential file (.pot)')
  parser.add_argument('-i','--input', help='in_struct.inp file name (if not specified, only pot file will be ploted)',required=False)

  parser.add_argument('-o','--output',type=str,help='Output file name for structure, by default use temp file and visualise it.', required=False)
  parser.add_argument('-O','--potential-output',type=str,help='Output file name for the potential. Default is a temp file, pass no argument to use <output_file>pot.<ext>', const=True, required=False, nargs='?')

  parser.add_argument('-f','--format',type=str,help='Output files format (see the allowed formats in ASE, default cif)', default='cif',required=False)
  parser.add_argument('-a','--visualise',help='Visualise the struct (default, if no output file is given)', action='store_true',required=False)
  parser.add_argument('-A','--visualise-potential',help='Visualise the potential', action='store_true',required=False)
  parser.add_argument('-s','--scale-radii',help='Ase visualisation atomic radius', type=float,required=False, default=0.5)
  parser.add_argument('-v','--vac',dest='vacuum_height', type=float,help='Visualisation-size of vacuum atoms in AA (default=10.0)',default=10.0,required=False)
  parser.add_argument('-b','--nbulk',type=int,help='Repetition of bulk unit (default=2)',default=2,required=False)


def run(args):
  from ...potentials.potentials import Potential          # NOQA
  from ...sprkkr.structure import structure_file_to_atoms # NOQA
  from ...ase.visualise import view                       # NOQA
  from tempfile import NamedTemporaryFile
  tmps = []

  def temp():
      file = NamedTemporaryFile()
      tmps.append(file)
      return file.name

  try:
      ciffile=args.output or temp()
      if args.potential_output is True:
          cf = ciffile.rsplit('.',1)
          ciffpotfile = '.'.join((cf[0] + '_pot', *cf[1:]))
      elif args.potential_output:
          ciffpotfile = args.potential_output
      else:
          ciffpotfile = temp()
      if not ciffile:
          ciffile = temp()

      outformat=args.format
      structure_filename=args.input

      # read the potential
      potential = Potential.from_file(args.potential_file)
      pot_atoms = potential.atoms

      def ase_view(atoms):
          view(atoms, scale_radii=args.scale_radii)

      if args.visualise_potential or (len(tmps)==2 and not structure_filename):
          ase_view(pot_atoms)
      pot_atoms.write(ciffpotfile, format = outformat)

      # From here on visualisation of in_structure.inp
      if structure_filename:
          # Extract data from structure file
          structure = structure_file_to_atoms(args.input, potential, n_bulk=args.nbulk, vacuum_height = args.vacuum_height)
          # Visualise the structure or write out to the file
          structure.write(ciffile, format = outformat)

          if args.visualise or (len(tmps)==2 and not args.visualise_potential):
              ase_view(structure)
  finally:
      for i in tmps:
          i.close()


if __name__ == "__main__":
    from ...common.tools import main
    main( globals() )
