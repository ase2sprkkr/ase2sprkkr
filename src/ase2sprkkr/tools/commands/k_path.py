#!/usr/bin/env python
"""
GUI for specifiing BSF k-points PATH
"""
from pathlib import Path
import sys
import argparse

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ...common.tools import main  # NOQA

help='GUI for specifiing BSF k-points PATH'

description="""\n Editor, where you can choose the points, along them is the k-path generated."""


def parser(parser):

    parser.add_argument('potential', help='Potential file.', type=str)
    parser.add_argument('-v', '--verbose', help='Verbose.', action='store_true', default=False)


def run(args):
    from ase2sprkkr import Potential
    from ase2sprkkr.gui.k_path import k_path_gui, k_path_to_string
    kpath = k_path_gui(Potential.from_file(args.potential).atoms, verbose=args.verbose)
    if kpath:
        print(k_path_to_string(kpath))
    else:
        print('No k-path selected')
        exit(-1)
