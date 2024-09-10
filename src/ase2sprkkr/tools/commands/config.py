#!/usr/bin/env python
"""
This command just run the tests
"""
from pathlib import Path
import sys

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help='Show the (current) user configuration.'
description='On the most modern unix/linux systems, it is in the file ~/.config/ase2sprkkr/__init__.py'


def parser(parser):
    parser.add_argument('-p', '--path', help='Just print the path to the file.', action='store_true')
    parser.add_argument('-i', '--info', help='Show the description of the configuration options.', action='store_true')
    parser.add_argument('-s', '--show', help='Print the configuration.', action='store_true')
    parser.add_argument('-e', '--edit', help='Edit the file using the editor in the $EDITOR environment variable.', action='store_true')
    parser.add_argument('-d', '--default', help='Put the default values into the file, if it not exists', action='store_true')
    parser.add_argument('-D', '--show-default', help='Show the default values', action='store_true')


def default_content(file):
    return f""" # ASE2SRPKKR configuration file
# -------------------------------

# Do not comment the following line.
from ase2sprkkr.config import config

# This file is pure python and it is executed when ase2sprkkr is imported.
# Place it into {file}

# This string is appended to the runned executables
# config.executables.suffix = ''

# Do you want to run the executables from a specific directory?
# config.executables.dir = ''

# Uncomment, if you don't want to run empty-spheres finding by default
# config.runing.empty_spheres = False

# Set to False if mpi should not be used. Or set to the number of processor,
# or just command line to run mpi programm, e.g.: [ '/usr/bin/mpiexec', '-n', '4' ]
# config.running.mpi = []

# You can change verbosity of the output setting to False or True
# config.runing.print_output = 'info'
"""


def run(args):
    import os
    from ...ase.register import user_preferences_file
    import subprocess
    from ...config import config
    file = user_preferences_file()

    run=True
    if args.default and not os.path.isfile(file):
        with open(file, 'w') as f:
            f.write(default_content(file))
        run=False
    if args.path:
        print(file)
        run=False
    if args.show_default:
        print(default_content(file))
        run=False
    if args.info:
        print(config._definition.description(verbose = 'all'))
        run=False
    if args.show:
        print(config.to_dict())
        run=False
    if args.edit:
        if 'EDITOR' not in os.environ:
            print("Please set the EDITOR environment variable.")
            exit(-1)
        subprocess.run([os.environ['EDITOR'], file])
        run=False
    if run:
        if os.path.isfile(file):
             print("# Content of the ASE2SPRKKR user configuration file")
             print("#--------------------------------------------------")
             print("")
             with open(file, 'r') as f:
                 print(f.read())
        else:
             print("# No custom user configuration found")


if __name__ == "__main__":
    main( globals() )
