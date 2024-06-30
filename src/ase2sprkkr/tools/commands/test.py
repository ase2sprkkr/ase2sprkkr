#!/usr/bin/env python
"""
This command just run the tests
"""
from pathlib import Path
import sys
import argparse
import os
import subprocess

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help='Run the tests of ASE2SPRKKR'
description='If something goes wrong, please send the output of this commands to developers'
unknowns = 'pytest_arguments'


def parser(parser):
    parser.add_argument('pytest_arguments', help='Arguments for pytest.', nargs=argparse.REMAINDER)
    parser.add_argument('--no-kkr', help='Do not run SPRKKR executables, just test the interface only.', action='store_false', default=True)


def run(args):
    if not args.no_kkr:
       os.environ['DO_NOT_RUN_SPRKKR'] = '1'
    a2s_path = os.path.join(root_path, 'ase2sprkkr')
    print(' '.join([sys.executable, '-m', 'pytest', '--doctest-modules', a2s_path] + args.pytest_arguments))
    subprocess.run([sys.executable, '-m', 'pytest', '--doctest-modules', a2s_path] + args.pytest_arguments)


if __name__ == "__main__":
    main( globals() )
