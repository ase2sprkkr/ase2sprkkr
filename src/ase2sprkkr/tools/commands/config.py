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
    parser.add_argument('-p', '--path', help='Just print the path to the configuration file.', action='store_true')
    parser.add_argument('-P', '--print', help='Just print the path configuration file.', action='store_true')
    parser.add_argument('-s', '--show', help='Print the configuration.', action='store_true')
    parser.add_argument('-i', '--info', help='Show the description of the configuration options.', action='store_true')
    parser.add_argument('-S', '--set', nargs=2,  help='Set the given configuration to the given value. Example: "ase2sprkkr config -S executables.suffix 8.6".', metavar=("NAME", "VALUE"))
    parser.add_argument('-e', '--edit', help='Edit the file using the editor in the $EDITOR environment variable.', action='store_true')
    parser.add_argument('-d', '--default', help='Put the default values into the file, if it not exists.', action='store_true')
    parser.add_argument('-D', '--overwrite-by-default', dest='default', help='Put the default values into the file. Owerwrite it if it exists.', action='store_const', const='overwrite')
    parser.add_argument('-o', '--show-default', help='Show the default values.', action='store_true')


def default_content(file):
    return f"""# ASE2SRPKKR configuration file
# -------------------------------

# Please, DO NOT comment the following line.
from ase2sprkkr.configuration import config

# This file is pure python and it is executed when ase2sprkkr is imported.
# Place it into the following path:
# {file}

# This string is appended to the runned executables
# config.executables.suffix = ''

# Do you want to run the executables from a specific directory?
# config.executables.dir = ''

# Uncomment, if you don't want to run empty-spheres finding by default
# config.running.empty_spheres = False

# Set to False if MPI should not be used. Or set to the number of processor,
# or just command line to run mpi programm, e.g.: [ '/usr/bin/mpiexec', '-n', '4' ]
# config.running.mpi = []

# You can change the verbosity of the output by setting the following to False or True
# config.running.print_output = 'info'

# Authentication token to Nomad. You can set it using ase2sprkkr nomad authenticate <username>
# config.nomad.token = None
"""


def run(args):
    import os
    import pyparsing
    from ...configuration import user_preferences_file, config
    import subprocess
    file = user_preferences_file()

    run=True

    if args.default:
        if args.default == 'overwrite' or not os.path.isfile(file):
            with open(file, 'w') as f:
                f.write(default_content(file))
            print(f"Configuration defaults have been written to file '{file}'.")
        else:
            print(f"Configuration file have already existed, so it have not been overwritten.")
        run=False
    if args.path:
        print(file)
        run=False
    if args.print:
        with open(file, 'r') as f:
           print(f.read())
        run=False
    if args.show_default:
        print(default_content(file))
        run=False
    if args.info:
        print(config._definition.description(verbose = 'all'))
        run=False
    if args.set:
        from ...common.grammar_types import Variant
        from ...common.warnings import DataValidityError
        try:
          try:
            val = Variant().parse(args.set[1])
            config.find(args.set[0]).set_permanent(val)
          except DataValidityError:
            #if parsing using mixed failed, just try store string
            config.find(args.set[0]).set_permanent(args.set[1])
          print(f"Configuration '{args.set[0]}' have been set to '{args.set[1]}'.")
        except KeyError:
          print(f"Unknown configuration name '{args.set[0]}'.")
        except pyparsing.ParseBaseException:
          print(f"Suspicious configuration value '{args.set[1]}'. If you are sure,"
                 "please edit the configuration file manually")
        run=False
    if args.show:
        import pprint
        pprint.pp(config.to_dict())
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
