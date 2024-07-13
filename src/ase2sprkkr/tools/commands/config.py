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

help='Show user configuration.'
description='On the most modern unix/linux systems, it is in the file ~/.config/ase2sprkkr/__init__.py'


def parser(parser):
    parser.add_argument('-p', '--path', help='Just print the path to the file.', action='store_true')
    parser.add_argument('-e', '--edit', help='Edit the file using the editor in the $EDITOR environment variable.', action='store_true')


def run(args):
    import os
    from ...ase.register import user_preferences_file
    import subprocess
    file = user_preferences_file()
    if args.path:
        print(file)
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
