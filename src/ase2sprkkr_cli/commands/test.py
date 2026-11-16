#!/usr/bin/env python
"""
This command just run the tests
"""

import sys


help = "Run the tests of ASE2SPRKKR."
description = "If something goes wrong, please send the output of this commands to developers to ase2sprkkr@ntc.zcu.cz"
unknowns = "pytest_arguments"


def parser(parser):
    # Pytest arguments must stay in parse_known_args()'s remainder.  Declaring
    # an argparse.REMAINDER positional here splits unknown options such as
    # ``-k EXPRESSION`` and loses their original ordering.
    parser.set_defaults(pytest_arguments=[])
    parser.add_argument(
        "--no-kkr",
        help="Do not run SPRKKR executables, just test the interface only.",
        action="store_false",
        default=True,
    )
    parser.add_argument("--pp", help="Pyparsing verbose stacktrace.", action="store_true")
    parser.add_argument("--fast", help="Exclude tests marked slow (useful locally; CI runs all tests).", action="store_true")
    smoke = parser.add_mutually_exclusive_group()
    smoke.add_argument("--smoke", dest="smoke", action="store_const", const=True,
                       default=None, help="Run only representative smoke tests.")
    smoke.add_argument("--no-smoke", dest="smoke", action="store_const", const=False,
                       help="Exclude smoke tests.")


def _pytest_selection(pytest_arguments, *, fast=False, smoke=None):
    """Compose CLI and forwarded pytest marker filters without overriding -m."""
    conditions = []
    if smoke is True:
        conditions.append("smoke")
    elif smoke is False:
        conditions.append("not smoke")
    if fast:
        conditions.append("not slow")
    if not conditions:
        return list(pytest_arguments)

    # Pytest's -m is a single-valued option: multiple occurrences silently
    # override one another. Intersect it with our filters instead.
    forwarded = []
    i = 0
    while i < len(pytest_arguments):
        arg = pytest_arguments[i]
        if arg in ("-m", "--markexpr"):
            if i + 1 == len(pytest_arguments):
                # Leave invalid options to pytest for its normal error handling.
                forwarded.append(arg)
                i += 1
                continue
            conditions.append(f"({pytest_arguments[i + 1]})")
            i += 2
            continue
        if arg.startswith(("-m=", "--markexpr=")):
            conditions.append(f"({arg.split('=', 1)[1]})")
        else:
            forwarded.append(arg)
        i += 1

    forwarded.extend(("-m", " and ".join(f"({value})" for value in conditions)))
    return forwarded


def run(args, global_args):
    from pathlib import Path
    import os
    import subprocess
    from pyparsing import ParserElement

    import ase2sprkkr, ase2sprkkr_cli

    if args.pp:
        ParserElement.verbose_stacktrace = True
    if not args.no_kkr:
        os.environ["DO_NOT_RUN_SPRKKR"] = "1"

    root_path = Path(ase2sprkkr.__file__).parent.parent
    packages = [ ase2sprkkr, ase2sprkkr_cli ]
    pytest_args = [
        sys.executable,
        "-m",
        "pytest",
        "--doctest-modules",
        "--rootdir",
        root_path,
        *( Path(i.__file__).parent for i in packages )
    ]
    pytest_args += _pytest_selection(args.pytest_arguments, fast=args.fast, smoke=args.smoke)
    out = subprocess.run(pytest_args)
    exit(out.returncode)

