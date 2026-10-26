#!/usr/bin/env python
"""
Show known structured result files, description of InputParameters and so on.
"""

from pathlib import Path
import sys
import argparse

if not __package__:
    __package__ = "ase2sprkkr.tools.commands"
sys.path.append(str(Path(__file__).resolve().parents[3]))


help = "Show capabilities of ase2sprkkr (known tasks, result file formats etc.)."
description = """\n Example of usage:
    #show all tasks
    ase2sprkkr info --task

    #show information about task SCF
    ase2sprkkr info --task=scf

    #show information about section ENERGY in task SCF
    ase2sprkkr info --task=scf.energy

    #find all sections/values with name ENERGY in all tasks
    ase2sprkkr info --task=.energy
"""


def parser(parser):
    parser.add_argument(
        "--task",
        "-t",
        help="Show informations about the known tasks. You can limit the output to a given task/section/option by an optional agrument [TASK]",
        nargs="?",
        const=True
    )
    parser.add_argument(
        "--output-file", "-o", help="Show information about known structured result files.", nargs="?", const=True
    )
    parser.add_argument(
        "--executables",
        action="store_true",
        help="Show known SPR-KKR executables and their default output readers.",
    )
    action = getattr(argparse, "BooleanOptionalAction", "store_true")
    parser.add_argument(
        "--verbose",
        "-v",
        help="Give a more verbose info. It is set automatically, if the output is limited.",
        action=action,
        default=None,
    )


def run(args, global_args):
    from ... import InputParameters, ResultFile
    from ...sprkkr.executables import executables  # NOQA
    from ...common.tools import main  # NOQA


    print("ASE2SPRKKR\n==========\nA tool to run SPRKKR calculations\n\n")

    def info_task(args, filter):
        registry = InputParameters.input_parameters_registry
        task_registry = InputParameters.task_registry

        if filter is True:
            path = False
            print(
                'The known SPRKKR tasks (kind of input parameters)\n'
                '---------------------------------------------------'
            )
            tasks = 'the known SPRKKR tasks'
            anyof = 'any of '
        else:
            if '.' in filter:
                filter, path = filter.split('.', 1)
            else:
                path = None


        entries = [
                (task.name, registry.definition_for_task(task), task.experimental)
                for task in task_registry.values(name_contains=filter.upper())
        ]

        if filter is not True:
            if not entries:
                print(f"There is no known SPRKKR task with '{filter.upper()}' in its name")
                return
            if filter == '':
                tasks = 'the known SPRKKR tasks'
                anyof = 'any of '
            elif len(entries) == 1 and entries[0][0].upper() == filter.upper():
                tasks = f"the task '{entries[0][0]}'"
                anyof = ''
            else:
                tasks = f"the tasks containing '{filter}'"
                anyof = 'any of '
            if path:
                title = f"The configuration sections/values '{path.upper()}' in {anyof}{tasks}"
            else:
                title = tasks.capitalize()
            print(f"{title}\n{'-' * len(title)}")

        # Keep the task identity distinct from its possibly shared schema.
        entries.sort(key=lambda entry: (entry[2], entry[0]))
        found = False
        experimental_started = False
        for task_name, definition, experimental in entries:
            members = None
            if path is not False:
                members = list(definition.create_object().get_members(path, is_option=None))
                if not members:
                    continue

            if experimental and not experimental_started:
                print("\n  Experimental tasks")
                print("  ~~~~~~~~~~~~~~~~~~")
                experimental_started = True

            if members is not None:
                for member in members:
                    found = True
                    member = member._definition
                    title = f"{member.item_type.capitalize()} {member.get_path()} of task {task_name}:"
                    print(title)
                    if args.verbose is not False:
                        print(member.description(prefix='  ', verbose=True))
                    else:
                        print('  ' + member.info())
            elif args.verbose:
                print(f"\nTask {task_name}:")
                print(definition.description(verbose=True))
            else:
                summary = definition.info()
                prefix = definition.display_name
                if summary.startswith(prefix + ' - '):
                    summary = summary[len(prefix) + 3:]
                elif summary == prefix:
                    summary = ''
                print(f"{task_name} - {summary}" if summary else task_name)

        if path and not found:
            print(f"There is no member named '{path}' in {anyof}{tasks}")

    def info_output_file(args, filter):
        print("Known SPR-KKR structured result files\n------------------------------------")
        if filter is not True:
            filter = filter.lower()
        for ext, i in ResultFile.definitions.items():
            if filter is not True and ((filter not in ext.lower()) or (filter not in i.name.lower())):
                continue

            if args.verbose:
                print("\n\n")
            print(f"{ext.upper()}: ", end="")

            if args.verbose or args.verbose is None and filter is not True:
                print(i.description(verbose=True))
            else:
                print(i.info())

    def info_executables():
        print("Known SPR-KKR executables\n-------------------------")
        definitions = tuple(executables.values())
        name_width = max(len("Executable"), *(len(item.name) for item in definitions))
        reader_width = max(
            len("Default reader"),
            *(len(item.default_output_reader) for item in definitions),
        )
        print(
            f"{'Executable':<{name_width}}  "
            f"{'Default reader':<{reader_width}}  Description"
        )
        for executable in definitions:
            print(
                f"{executable.name:<{name_width}}  "
                f"{executable.default_output_reader:<{reader_width}}  "
                f"{executable.description}"
            )

    printed = False
    if args.task is not None:
        filter = args.task.lower() if isinstance(args.task, str) else args.task
        info_task(args, filter)
        print("\n")
        printed = True
    if args.output_file is not None:
        filter = args.output_file.lower() if isinstance(args.output_file, str) else args.output_file
        info_output_file(args, filter)
        print("\n")
        printed = True
    if args.executables:
        info_executables()
        print("\n")
        printed = True

    if not printed:
        info_task(args, True)
        print("\n")
        info_output_file(args, True)
        print("\n")
