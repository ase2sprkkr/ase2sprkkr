#!/usr/bin/env python
"""
This command just show the path to the examples
"""
from pathlib import Path
import sys
import os

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help="Show path to the examples, or example itself, or copy a given example."
description = "See also the 'shell -e' subcommand for interactivelly running the example"

def parser(parser):
    parser.add_argument('example', type=int, nargs='?', help='The number of the example to print. If ommited')
    parser.add_argument('-c', '--copy', help='Copy the example to a given dir', type=str)
    parser.add_argument('-p', '--path', help='Print a path to the example(s)', action='store_true')

def run(args):
    from ase2sprkkr.gui import examples

    if not args.example:
        if not args.path:
            print("Examples dir: ", end='')
        print(examples.examples_dir())
        if args.path:
            return
        print("")

        exs = examples.list_of_examples()
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
    if show:
        print(example.source())


if __name__ == "__main__":
    main( globals() )
