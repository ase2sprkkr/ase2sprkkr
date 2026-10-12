#!/usr/bin/env python
"""Tool for ASE2SPRKKR user configuration."""

from pathlib import Path
import sys

if not __package__:
    __package__ = "ase2sprkkr.tools.commands"

root_path = str(Path(__file__).resolve().parents[3])
sys.path.append(root_path)

from ...common.tools import main  # NOQA

help = "Show or edit the current user configuration."
description = "The user configuration is stored in ~/.config/ase2sprkkr/config.yaml on typical Unix/Linux systems."


def parser(parser):
    parser.add_argument("-p", "--path", help="Just print the path to the configuration file.", action="store_true")
    parser.add_argument("-P", "--print", help="Print the configuration file.", action="store_true")
    parser.add_argument("-s", "--show", help="Print the resolved configuration.", action="store_true")
    parser.add_argument("-i", "--info", help="Show the description of the configuration options.", action="store_true")
    parser.add_argument(
        "-S",
        "--set",
        nargs=2,
        help='Set the given configuration to the given value. Example: "ase2sprkkr config -S default_profile local".',
        metavar=("NAME", "VALUE"),
    )
    parser.add_argument(
        "-e", "--edit", help="Edit the YAML file using the editor in the $EDITOR environment variable.", action="store_true"
    )
    parser.add_argument(
        "-d", "--default", help="Create a commented YAML configuration file if it does not exist.", action="store_true"
    )
    parser.add_argument(
        "-D",
        "--overwrite-by-default",
        dest="default",
        help="Write the commented default YAML configuration, overwriting an existing file.",
        action="store_const",
        const="overwrite",
    )
    parser.add_argument("-o", "--show-default", help="Show the default YAML configuration.", action="store_true")


def default_content(file):
    return f"""# ASE2SPRKKR configuration file
# -----------------------------
# Location: {file}
#
# Only values you want to override need to be present. The examples below are
# commented out and therefore do not change the defaults.

# Explicitly select one of the named profiles below.
# default_profile: local

# profiles:
#   local:
#     executable_dir: /opt/sprkkr/bin
#     executable_suffix: ""
#     # MPI can be false, true, an integer process count, "auto", or argv.
#     mpi: auto
#
# The reserved profile "ase" is read from ASE configuration and is not stored
# here. The reserved "legacy" profile is managed through deprecated options and
# may appear here after migrating an old configuration.

# running:
#   empty_spheres: auto
#   print_output: info
#   mpi_warning: true

# nomad:
#   token: your-token
"""


def run(args, global_args):
    import os
    import pyparsing
    from ...configuration import user_preferences_file, config
    import subprocess

    file = user_preferences_file()

    run = True

    if args.default:
        if args.default == "overwrite" or not os.path.isfile(file):
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            with open(file, "w", encoding="utf8") as f:
                f.write(default_content(file))
            print(f"Configuration defaults have been written to file '{file}'.")
        else:
            print("Configuration file already exists, so it has not been overwritten.")
        run = False
    if args.path:
        print(file)
        run = False
    if args.print:
        with open(file, "r", encoding="utf8") as f:
            print(f.read())
        run = False
    if args.show_default:
        print(default_content(file))
        run = False
    if args.info:
        print(config._definition.description(verbose="all"))
        run = False
    if args.set:
        from ...common.grammar_types import Variant
        from ...common.warnings import DataValidityError

        try:
            try:
                val = Variant().parse(args.set[1])
                config.find(args.set[0]).set_permanent(val)
            except DataValidityError:
                # If parsing using mixed failed, just try storing a string.
                config.find(args.set[0]).set_permanent(args.set[1])
            print(f"Configuration '{args.set[0]}' has been set to '{args.set[1]}'.")
        except KeyError:
            print(f"Unknown configuration name '{args.set[0]}'.")
        except pyparsing.ParseBaseException:
            print(
                f"Suspicious configuration value '{args.set[1]}'. If you are sure, "
                "please edit the configuration file manually."
            )
        run = False
    if args.show:
        import pprint

        pprint.pp(config.to_dict())
        run = False
    if args.edit:
        if "EDITOR" not in os.environ:
            print("Please set the EDITOR environment variable.")
            raise SystemExit(-1)
        if not os.path.isfile(file):
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            Path(file).write_text(default_content(file), encoding="utf8")
        subprocess.run([os.environ["EDITOR"], file])
        run = False

    if run:
        if os.path.isfile(file):
            print("# Content of the ASE2SPRKKR user configuration file")
            print("#--------------------------------------------------")
            print("")
            with open(file, "r", encoding="utf8") as f:
                print(f.read())
        else:
            print("# No custom user configuration found")


if __name__ == "__main__":
    main(globals())
