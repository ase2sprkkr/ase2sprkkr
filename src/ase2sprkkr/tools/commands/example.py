#!/usr/bin/env python
"""
This command just show the path to the examples
"""

from pathlib import Path
import sys

if not __package__:
    __package__ = "ase2sprkkr.tools.commands"

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help = "Show path to the examples, or example itself, or copy a given example."
description = "See also the 'shell -e' subcommand for running examples interactivelly"
aliases = ["examples"]


def parser(parser):
    parser.add_argument("example", type=int, nargs="?", help="The number of the example to print. If ommited")
    parser.add_argument("-c", "--copy", help="Copy the example to a given dir", type=str)
    parser.add_argument("-p", "--path", help="Print a path to the example(s)", action="store_true")
    parser.add_argument("-s", "--script", help="Print a path to the example main script", action="store_true")
    parser.add_argument("-g", "--regex", help="Find an example according to the given regex", type=str)


def run(args, global_args):
    from ase2sprkkr.gui import examples

    if not args.example:
        if args.script:
            raise ValueError("The --script argument requires to specify the example number.")
        if not args.path:
            print("Examples dir: ", end="")
        print(examples.examples_dir())
        if args.path:
            return
        print("")

        exs = examples.list_of_examples(regex=args.regex)
        print(f"{'NAME':<30} {'SCRIPT':<20} DESCRIPTION")
        for e in exs:
            print(f"{e.name:<30} {e.main_script.name:<20} {e.short_docstring or ''}")
        return

    example = examples.Example.by_number(args.example)

    show = True
    if args.copy:
        example.copy(args.copy)
        show = False
    if args.path:
        print(example.dir)
        show = False
    if args.script:
        print(example.main_script)
        show = False
    if show:
        print(example.source())


if __name__ == "__main__":
    main(globals())
