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
    parser.add_argument("--fast", help="Run only tests not marked as slow.", action="store_true")


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
    if args.fast:
        pytest_args += ["-m", "not slow"]
    pytest_args += args.pytest_arguments
    out = subprocess.run(pytest_args)
    exit(out.returncode)

