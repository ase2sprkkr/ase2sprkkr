#!/usr/bin/env python
"""
GUI for specifiing BSF k-points PATH
"""

from pathlib import Path
import sys

if not __package__:
    __package__ = "ase2sprkkr.tools.commands"
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ...common.tools import main  # NOQA

help = "GUI for specifiing BSF k-points PATH (for E-K mode)"

description = """\n Editor, where you can choose the points, along them is the k-path generated. The output can be added to the input parameters, e.g.

from numpy import array
calculator.input_parameters.set(
{'NKDIR': 1, 'NK': 146, 'KA': array([[ 0.33333333, -0.57735027, -0.02551552]]), 'KE': array([[ 6.66666667e-01, -8.21078625e-17, -2.55155182e-02]])}
)
"""


def parser(parser):

    parser.add_argument("potential", help="Potential file.", type=str)
    parser.add_argument(
        "-v", "--verbose", help="Verbose.", action="store_true", default=False
    )
    parser.add_argument(
        "-i", "--input-file", help="Print the path as part of the SPRKKR input file."
    )


def run(args, global_args):
    from ase2sprkkr import Potential
    from ase2sprkkr.gui.k_path import k_path_gui, k_path_to_string

    kpath = k_path_gui(Potential.from_file(args.potential).atoms, verbose=args.verbose)
    if kpath:
        if args.input_file:
            print(k_path_to_string(kpath))
        else:
            print(kpath)
    else:
        print("No k-path selected")
        exit(-1)
