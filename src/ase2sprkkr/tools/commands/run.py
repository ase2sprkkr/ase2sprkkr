#!/usr/bin/env python
"""
This is a sctipt to visualise in_struct.inp files. Run it to see the doc.
"""
from pathlib import Path
import sys
import argparse

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ... import InputParameters, SPRKKR, Potential  # NOQA
from ...common.tools import parse_named_option, main  # NOQA


help='Run SPRKKR calculation.'
description='You can pass the input parameters using the command line. However, due to limitation of argparse, pass them directly after the potential file.'


def parser(parser):
    parser.add_argument('pot', help='SPR-KKR potential file')
    parser.add_argument('-t','--task', dest='task', choices=InputParameters.definitions.keys(), type = str.lower, help='Task to compute')
    parser.add_argument('--print-output', '-O', help='Print the output of SPRKKR', action=argparse.BooleanOptionalAction)
    parser.add_argument('--output-file', '-o', help='Output file', type=str)
    parser.add_argument('--input-file', '-i', help='Input (parameters) file', type=str)
    parser.add_argument('--empty-spheres', '-e', help='Try to add empty spheres (by default, they are added if the potential is not converged and no vacuum atom is present', action=argparse.BooleanOptionalAction)
    parser.add_argument('options', nargs='*', help='Input parameters in the form <name>=<value>', type=parse_named_option)


def run(args):
    pot = Potential.from_file(args.pot)
    calc = SPRKKR(potential=pot)
    if args.task:
        calc.input_parameters = args.task
    if args.input_file:
        calc.input_parameters = args.input_file
    kwargs = {}
    for i in 'empty_spheres', 'print_output', 'output_file':
        o=getattr(args, i)
        if o is not None:
            kwargs[i]=o
    calc.calculate(options=dict(args.options), **kwargs)


if __name__ == "__main__":
    main( globals() )