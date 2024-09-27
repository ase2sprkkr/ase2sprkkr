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

from ...common.tools import parse_named_option, main  # NOQA
from ... import InputParameters  # NOQA


help='Run SPRKKR calculation.'
description='You can pass the input parameters using the command line. However, due to limitation of argparse, pass them directly after the potential file.'


def parser(parser):
    bool_arg = getattr(argparse, 'BooleanOptionalAction', 'store_true')
    parser.add_argument('pot', help='SPR-KKR potential file or input file. The type is determined by its content.')

    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-t','--task', dest='task', choices=InputParameters.definition_modules.keys(), type = str.lower, help='Task to compute.')
    grp.add_argument('--input-file', '-i', help='Run ase2sprkkr -i inputfile potfile, if you want to run task determined by the give input file with a given potential.', type=str)

    parser.add_argument('--print-output', '-O', help='Print the output of SPRKKR.', default='None', action=bool_arg)
    parser.add_argument('--output-file', '-o', help='Output file', type=str)
    parser.add_argument('--empty-spheres', '-e', help='Try to add empty spheres (by default, they are added if the potential is not converged and no vacuum atom is present.', action=bool_arg)
    parser.add_argument('options', nargs='*', help='Input parameters in the form <name>=<value>.', type=parse_named_option)


def run(args):
    from ... import SPRKKR, Potential  # NOQA
    import pyparsing
    import sys

    calc = SPRKKR()
    try:
        potential=Potential.from_file(args.pot)
    except pyparsing.ParseException:
        try:
            calc.input_parameters = ip = InputParameters.from_file(args.pot)
            potential = ip.CONTROL.POTFIL()
        except pyparsing.ParseException:
            print(f"File {args.pot} is neither potential nor input_parameters file")
            sys.exit(-1)

    if args.input_file:
        calc.input_parameters = args.input_file
    elif args.task:
        calc.input_parameters = args.task
    kwargs = {}
    for i in 'empty_spheres', 'print_output', 'output_file':
        o=getattr(args, i)
        if o is not None:
            kwargs[i]=o

    calc.calculate(potential=potential,
                   input_parameters=ip,
                   options=dict(args.options), **kwargs)


if __name__ == "__main__":
    main( globals() )
