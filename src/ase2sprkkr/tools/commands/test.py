#!/usr/bin/env python
"""
This command just run the tests
"""
from pathlib import Path
import sys
import argparse

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help='Run the tests of ASE2SPRKKR.'
description='If something goes wrong, please send the output of this commands to developers to ase2sprkkr@ntc.zcu.cz'
unknowns = 'pytest_arguments'


def parser(parser):
    parser.add_argument('pytest_arguments', help='Arguments for pytest.', nargs=argparse.REMAINDER)
    parser.add_argument('--no-kkr', help='Do not run SPRKKR executables, just test the interface only.', action='store_false', default=True)
    parser.add_argument('--pp', help='Pyparsing verbose stacktrace.', action='store_true')
    parser.add_argument('--fast', help='Run only tests not marked as slow.', action='store_true')


def run(args, global_args):
    import pytest
    import os
    import subprocess
    from pyparsing import ParserElement
    if args.pp:
        ParserElement.verbose_stacktrace=True
    if not args.no_kkr:
       os.environ['DO_NOT_RUN_SPRKKR'] = '1'
    a2s_path = os.path.join(root_path, 'ase2sprkkr')
    if pytest.version_tuple[0] <= 6:
        print("Pytest version >= 6.0 is required. Please install it with pip install --upgrade 'pytest>=6'")
        exit(-1)
    pytest_args = [sys.executable, '-m', 'pytest', '--doctest-modules', '--rootdir', root_path, a2s_path]
    if args.fast:
        pytest_args += ['-m', 'not slow']
    pytest_args += args.pytest_arguments
    out = subprocess.run(pytest_args)
    exit(out.returncode)


if __name__ == "__main__":
    main( globals() )
