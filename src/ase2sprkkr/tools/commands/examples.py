#!/usr/bin/env python
"""
This command just run the tests
"""
from pathlib import Path
import sys
import os

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help='Show path to the examples.'
description=''


def parser(parser):
    pass


def run(args):
    import ase2sprkkr.examples as e
    print(os.path.dirname(e.__file__))


if __name__ == "__main__":
    main( globals() )
