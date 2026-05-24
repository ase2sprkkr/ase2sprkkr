#!/usr/bin/env python
"""
Opening the ASE2SPRKKR output files and/or scripts in the shell or Jupyter notebook
"""
from pathlib import Path
import sys
import argparse

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ...common.tools import parse_tuple_function, parse_named_option, \
        append_id_to_filename, parse_inches, main # NOQA
from ...common.lazy_string import LazyString      # NOQA


@LazyString
def description():
    from ...output_files.output_files import OutputFile
    out = 'The type of the file is guessed from the content of the file and from the extension. The currently supported files are: \n'
    defs = OutputFile.definitions.items()
    out += '\n'.join(map(lambda x: f"    {x[0].upper()}: {x[1].info()}", defs))
    return out


help='Open SPR-KKR output files in a Python shell or Jupyter notebook'

def parser(parser):
    parser.add_argument('-o', '--output', help='SPR-KKR output file name (see the supported files above).', action='append', default=[])
    parser.add_argument('-j', '--jupyter', help='Open Jupyter notebook.', action="store_true")
    parser.add_argument('-D', '--pdb', help='Run the command inside pdb', action="store_true")
    parser.add_argument('-e', '--example', help='Open and run the example given by number. The example files will be copied to the directory given by -d or to a temp directory.', action="store", type=int)
    parser.add_argument('-d', '--directory', help='Directory, where to copy the example. Use "." for the current directory.')
    parser.add_argument('-s', '--save', help='Save the script/Jupyter notebook.', action="store", type=str)
    parser.add_argument('-r', '--run', help='Execute the notebook before saving and opening Jupyter notebook.', action="store", type=str)

def run(args, global_args):
  from ...gui.shell import Python, JupyterLab, Pdb
  import tempfile
  import contextlib

  if args.jupyter:
    shell = JupyterLab
    kwargs = { 'run': args.run }
  elif args.pdb:
    shell = Pdb
  else:
    shell = Python
  kwargs = {}
  shell = shell(**kwargs)

  full_path = True
  if args.directory:
      dire = contextlib.nullcontext(args.directory)
  elif args.example:
      dire = tempfile.TemporaryDirectory()
  else:
      full_path = False
      dire = contextlib.nullcontext(None)

  for i in args.output:
      if full_path:
          i = Path(file_path).resolve(i)
      shell.add_output_file(i)

  with dire as dr:
      if dr is not None:
          shell.change_working_dir(dr)

      if args.example:
          shell.run_example(args.example, dr)

      if args.save:
          shell.save(args.save)

      shell.open()

if __name__ == "__main__":
    main( globals() )
