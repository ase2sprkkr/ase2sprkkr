#!/usr/bin/env python
"""
This is a script to run calculation from commandline. Run it to see the doc.
"""

import argparse

from ase2sprkkr_cli.errors import CLIUsageError




help = "Run SPRKKR calculation."
description = "You can pass the input parameters using the command line. However, due to limitation of argparse, pass them directly after the potential file."


def parser(parser):
    bool_arg = getattr(argparse, "BooleanOptionalAction", "store_true")
    parser.add_argument("pot", help="SPR-KKR potential file or input file. The type is determined by its content.")

    grp = parser.add_mutually_exclusive_group()
    grp.add_argument(
        "-t",
        "--task",
        dest="task",
        type=str,
        help="Task to compute (see 'ase2sprkkr info --task').",
    )
    grp.add_argument(
        "--input-file",
        "-i",
        help="Run ase2sprkkr -i inputfile potfile, if you want to run task determined by the give input file with a given potential.",
        type=str,
    )

    parser.add_argument("--print-output", "-O", help="Print the output of SPRKKR.", default="None", action=bool_arg)
    parser.add_argument("--output-file", "-o", help="Output file", type=str)
    parser.add_argument(
        "--empty-spheres",
        "-e",
        help="Try to add empty spheres (by default, they are added if the potential is not converged and no vacuum atom is present.",
        action=bool_arg,
    )
    parser.add_argument(
        "options", nargs="*", help="Input parameters in the form <name>=<value>.", type=str
    )


def run(args, global_args):
    from ase2sprkkr import InputParameters
    from ase2sprkkr.common.module_registry import RegistryEntryNotFound
    from ase2sprkkr.common.tools import parse_named_option
    import pyparsing

    if args.task is not None:
        try:
            args.task = InputParameters.task_registry[args.task].name.lower()
        except RegistryEntryNotFound:
            raise CLIUsageError(f"Unknown task {args.task!r}") from None
    try:
        options = dict(parse_named_option(option) for option in args.options)
    except (ValueError, pyparsing.ParseBaseException) as exc:
        raise CLIUsageError(f"Invalid input parameter: {exc}") from None
    from ase2sprkkr import SPRKKR, Potential  # NOQA
    import sys

    calc = SPRKKR()
    ip = None
    try:
        potential = Potential.from_file(args.pot)
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
    for i in "empty_spheres", "print_output", "output_file":
        o = getattr(args, i)
        if o is not None:
            kwargs[i] = o

    calc.calculate(potential=potential, input_parameters=ip, options=options, **kwargs)

