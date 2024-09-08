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
    parser.add_argument('-e', '--edit', help='Edit the file using the editor in the $EDITOR environment variable.', action='store_true')
    parser.add_argument('-d', '--default', help='Put the default values into the file, if it not exists', action='store_true')
    parser.add_argument('-D', '--show-default', help='Show the default values', action='store_true')


def default_content(file):
    return f""" # ASE2SRPKKR configuration file
# -------------------------------

# Do not comment the following line.
from ase2sprkkr import config

# This file is pure python and it is executed when ase2sprkkr is imported.
# Place it into {file}

# This string is appended to the runned executables
# config.sprkkr_executable_suffix = ''

# Uncomment, if you don't want to run empty-spheres finding by default
# config.empty_spheres = False
"""


def run(args):
    import os
    from ...ase.register import user_preferences_file
    import subprocess
    file = user_preferences_file()

    if args.default and not os.path.isfile(file):
        with open(file, 'w') as f:
            f.write(default_content(file))
    if args.path:
        print(file)
    elif args.show_default:
        print(default_content(file))
    elif args.edit:
        if 'EDITOR' not in os.environ:
            print("Please set the EDITOR environment variable.")
            exit(-1)
        subprocess.run([os.environ['EDITOR'], file])
    else:
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
