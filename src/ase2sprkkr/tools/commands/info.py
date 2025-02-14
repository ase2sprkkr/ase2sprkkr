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
description="""\n Example of usage:
    #show all tasks
    ase2sprkkr info --task

    #show information about task SCF
    ase2sprkkr info --task=scf

    #show information about section ENERGY in task SCF
    ase2sprkkr info --task=scf.energy

    #find all sections/values with name ENERGY in all tasks
    ase2sprkkr info --task=.energy
"""


""" This options will be printed, if no argument is given """
default_options = [
  'task',
  'output_file',
]

""" These are the all options """
possible_options = default_options


def parser(parser):
    parser.add_argument('--task', '-t', help='Show informations about the known tasks. You can limit the output to a given task/section/option by an optional agrument [TASK]', nargs='?', const=True)
    parser.add_argument('--output-file', '-o', help='Show informations about the known output files. You can ', nargs='?', const=True)
    action = getattr(argparse, 'BooleanOptionalAction', 'store_true')
    parser.add_argument('--verbose', '-v', help='Give a more verbose info. It is set automatically, if the output is limited.', action=action, default=None)


def run(args):

    print("ASE2SPRKKR\n"
          "==========\n"
          "A tool to run SPRKKR calculations\n\n")

    def info_task(args, filter):
        input_parameters = InputParameters.definitions.values()
        if filter is True:
            input_parameters = list(input_parameters)
            path = None
            print("The known SPRKKR tasks (kind of input parameters)\n"
                  "---------------------------------------------------")
            tasks = 'the known SPRKKR tasks'
            anyof = 'any of '
        else:
            if '.' in filter:
                filter,path = filter.split('.', 1)
            else:
                path = None
            input_parameters = [ i for i in input_parameters if filter in i.name.lower() ]
            if len(input_parameters) == 0:
                print("There is no known SPRKKR task with '{filter.upper()}' in its name")
                return
            if filter == '':
                tasks = 'the known SPRKKR tasks'
                anyof = 'any of '
            elif len(input_parameters) == 1 and \
               input_parameters[0].name.upper() == filter.upper():
                   tasks = f"the task '{input_parameters[0].name.upper()}'"
                   anyof = ''
            else:
                tasks = f"the tasks containing '{filter}'"
                anyof = 'any of '
            if path:
                x=f"The configuration sections/values '{path.upper()}' in {anyof}{tasks}"
                print(f"{x}\n{'-'*len(x)}")
            else:
                x=f"{tasks.capitalize()}'"
                print(f"{x}\n{'-'*len(x)}")

        found = False

        for i in input_parameters:
            def print_ip(i, prefix='', add=None):
                if args.verbose or args.verbose is None and filter is not True:
                    print("\n")
                    if add:
                        print(add)
                    print(i.description(prefix=prefix, verbose=True))
                else:
                    if add:
                        print(add)
                    print(prefix + i.info())

            if path is not None:
                ii = list(i.create_object().get_members(path, is_option=False))
                if not ii:
                    continue
                found = True
                for val in ii:
                    val = val._definition
                    add = f"{val.item_type.capitalize()} {val.get_path()} of task {i.name.upper()}:"
                    print_ip(val, '  ', add)
            else:
                print_ip(i, prefix='')

        if path and not found:
            out = f"There is no member named '{path}' in {anyof}{tasks}"
            print(out)

    def info_output_file(args, filter):
        print("Known SPRKKR output files\n"
              "--------------------")
        if filter is not True:
            filter = filter.lower()
        for ext, i in OutputFile.definitions.items():
            i = i.definition
            if filter is not True and (
                (filter not in ext.lower()) or
                (filter not in i.name.lower())
              ):
                continue

            if args.verbose:
                print("\n\n")
            print(f"{ext.upper()}: ",end='')

            if args.verbose or args.verbose is None and filter is not True:
                print(i.description(verbose=True))
            else:
                print(i.info())

    g = locals()

    def print_option(i, filter):
        if filter is not True:
            filter = filter.lower()
        g["info_" + i](args, filter)
        print("\n")

    printed=False
    for i in possible_options:
        if getattr(args, i) is not None:
            print_option(i, getattr(args, i))
            printed=True
    if not printed:
        for i in default_options:
            print_option(i, True)
