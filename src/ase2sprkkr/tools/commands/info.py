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

from ... import InputParameters, OutputFile  # NOQA
from ...common.tools import main  # NOQA


help='Show capabilities of ase2sprkkr (known tasks, output files formats etc.).'
description='You can pass the input parameters using the command line. However, due to limitation of argparse, pass them directly after the potential file.'

always_options = [
  'tasks',
  'output_files',
]

possible_options = [] + always_options


class TestsArgAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        all_tests = possible_options
        default_tests = always_options

        if not values:
            setattr(namespace, self.dest, default_tests)
            return

        # If no argument is specified, the default gets passed as a
        # string 'default' instead of as a list ['default']. Probably
        # a bug in argparse. The below gives us a list.
        if not isinstance(values, list):
            values = [values]

        tests = set(values)

        # If 'all', is found, replace it with the tests it represents.
        # For reasons of compatibility, 'all' does not actually include
        # one of the tests (let's call it 'e'). So we can't just do
        # tests = all_tests.
        try:
            tests.remove('all')
            tests.update(set(all_tests))
        except KeyError:
            pass

        setattr(namespace, self.dest, sorted(list(tests)))


def parser(parser):
    parser.add_argument('about',
                        nargs='*',
                        help='About what do you want an information?',
                        metavar=f'{{{",".join(possible_options)}}}',
                        type=str,
                        default=[],
                        action=TestsArgAction
                        )
    parser.add_argument('--verbose', '-v', help='Give a more verbose info.', action='store_true')


def run(args):

    print("ASE2SPRKKR\n"
          "==========\n"
          "A tool to run SPRKKR calculations\n\n")

    g = globals()
    for i in always_options:
        if not args.about or i in args.about:
            g["info_" + i](args)
        print("\n")


def info_tasks(args):
    print("Possible task (kind of input parameters)\n"
          "----------------------------------------")
    for i in InputParameters.definitions.values():
        if args.verbose:
            print("\n\n")
            print(i.description(verbose=True))
        else:
            print(i.info())


def info_output_files(args):
    print("Known output files\n"
          "--------------------")
    for ext, i in OutputFile.definitions.items():
        if args.verbose:
            print("\n\n")
        print(f"{ext.upper()}: ",end='')
        i = i.definition
        if args.verbose:
            print(i.description(verbose=True))
        else:
            print(i.info())
