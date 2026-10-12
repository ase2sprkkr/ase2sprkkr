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
    parser.add_argument("-p", "--path", help="Print the path to config.yaml.", action="store_true")
    parser.add_argument("-P", "--print", help="Print config.yaml.", action="store_true")
    parser.add_argument("-s", "--show", help="Print the resolved configuration.", action="store_true")
    parser.add_argument("-i", "--info", help="Show the description of the configuration options.", action="store_true")
    parser.add_argument(
        "-S",
        "--set",
        nargs=2,
        help='Set a configuration value and persist it to YAML. Example: "ase2sprkkr config -S default_profile local".',
        metavar=("NAME", "VALUE"),
    )
    parser.add_argument(
        "-e", "--edit", help="Edit config.yaml using the editor in $EDITOR.", action="store_true"
    )
    parser.add_argument(
        "-d", "--default", help="Create a commented example config.yaml if it does not exist.", action="store_true"
    )
    parser.add_argument(
        "-D",
        "--overwrite-by-default",
        dest="default",
        help="Write the commented example config.yaml, overwriting an existing file.",
        action="store_const",
        const="overwrite",
    )
    parser.add_argument("-o", "--show-default", help="Show the commented example YAML configuration.", action="store_true")
    parser.add_argument(
        "-m",
        "--migrate",
        action="store_true",
        help="Migrate the deprecated Python user configuration to config.yaml.",
    )
    parser.add_argument(
        "--init-script",
        action="store_true",
        help="Print the path to the deprecated Python __init__.py configuration.",
    )
    parser.add_argument(
        "--print-init-script",
        action="store_true",
        help="Print the deprecated Python __init__.py configuration.",
    )
    parser.add_argument(
        "--edit-init-script",
        action="store_true",
        help="Edit the deprecated Python __init__.py configuration using $EDITOR.",
    )


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


def _edit(file, *, create=None):
    import os
    import subprocess

    editor = os.environ.get("EDITOR")
    if not editor:
        print("Please set the EDITOR environment variable.")
        raise SystemExit(-1)
    if create is not None and not Path(file).is_file():
        Path(file).parent.mkdir(parents=True, exist_ok=True)
        Path(file).write_text(create, encoding="utf8")
    subprocess.run([editor, file])


def run(args, global_args):
    import os
    import pyparsing
    from ...configuration import (
        config,
        legacy_user_preferences_file,
        user_preferences_file,
    )

    file = user_preferences_file()
    legacy_file = legacy_user_preferences_file()
    run = True

    if args.migrate:
        if os.path.isfile(file):
            print(f"YAML configuration already exists: '{file}'. Nothing was migrated.")
        elif not os.path.isfile(legacy_file):
            print(f"No legacy configuration found at '{legacy_file}'.")
        else:
            config.save_user_preferences()
            print(f"Migrated '{legacy_file}' to '{file}'. The legacy file was left unchanged.")
        run = False
    if args.default:
        if args.default == "overwrite" or not os.path.isfile(file):
            Path(file).parent.mkdir(parents=True, exist_ok=True)
            Path(file).write_text(default_content(file), encoding="utf8")
            print(f"Example configuration has been written to '{file}'.")
        else:
            print("Configuration file already exists, so it has not been overwritten.")
        run = False
    if args.path:
        print(file)
        run = False
    if args.print:
        if os.path.isfile(file):
            print(Path(file).read_text(encoding="utf8"), end="")
        else:
            print(f"No YAML configuration found at '{file}'.")
        run = False
    if args.show_default:
        print(default_content(file), end="")
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
        _edit(file, create=default_content(file))
        run = False
    if args.init_script:
        print(legacy_file)
        run = False
    if args.print_init_script:
        if os.path.isfile(legacy_file):
            print(Path(legacy_file).read_text(encoding="utf8"), end="")
        else:
            print(f"No legacy init script found at '{legacy_file}'.")
        run = False
    if args.edit_init_script:
        if not os.path.isfile(legacy_file):
            print(f"No legacy init script found at '{legacy_file}'.")
        else:
            _edit(legacy_file)
        run = False

    if run:
        if os.path.isfile(file):
            print("# Content of the ASE2SPRKKR YAML configuration")
            print("#---------------------------------------------")
            print()
            print(Path(file).read_text(encoding="utf8"), end="")
        elif os.path.isfile(legacy_file):
            print("No YAML configuration exists.")
            print()
            print(f"Legacy Python configuration found at:\n    {legacy_file}")
            print()
            print("Run 'ase2sprkkr config --migrate' to migrate it to config.yaml.")
        else:
            print("# No custom user configuration found")


if __name__ == "__main__":
    main(globals())
