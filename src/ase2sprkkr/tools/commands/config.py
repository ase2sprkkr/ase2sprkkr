#!/usr/bin/env python
"""Tool for ASE2SPRKKR user configuration."""

from pathlib import Path
from collections.abc import Mapping
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
        help=(
            "Set a configuration value and persist it to YAML. "
            "Examples: '-S default_profile local', "
            "'-S profiles.cluster.mpi 8', or "
            "'-S profiles.cluster \"{mpi: 8, executable_dir: /opt/sprkkr}\"'."
        ),
        metavar=("NAME", "VALUE"),
    )
    parser.add_argument(
        "-U",
        "--unset",
        metavar="NAME",
        help=(
            "Unset a configuration value and persist the change. "
            "For a whole named execution profile use, for example, "
            "'--unset profiles.cluster'."
        ),
    )
    parser.add_argument(
        "--list-profiles",
        action="store_true",
        help="List available execution profiles and mark the explicit default.",
    )
    parser.add_argument(
        "--executables",
        action="store_true",
        help="List known SPR-KKR executables and their default output readers.",
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
        help=(
            "Migrate declarative settings from the Python configuration hook to "
            "config.yaml, convert legacy execution settings to a named profile, "
            "and remove migrated top-level assignments from the hook."
        ),
    )
    cleanup = parser.add_mutually_exclusive_group()
    cleanup.add_argument(
        "--keep-init",
        action="store_true",
        help=(
            "After migration, keep the Python hook unchanged instead of removing "
            "migrated assignments."
        ),
    )
    cleanup.add_argument(
        "--clean-init",
        action="store_true",
        help=(
            "Remove migrated top-level ASE2SPRKKR assignments from the Python hook "
            "after migration (the default)."
        ),
    )
    cleanup.add_argument(
        "--delete-init",
        action="store_true",
        help="After migration, delete the Python configuration hook completely.",
    )
    parser.add_argument(
        "--init-script",
        action="store_true",
        help="Print the path to the Python user-configuration hook.",
    )
    parser.add_argument(
        "--print-init-script",
        action="store_true",
        help="Print the Python user-configuration hook.",
    )
    parser.add_argument(
        "--edit-init-script",
        action="store_true",
        help="Edit the Python user-configuration hook using $EDITOR.",
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
#     # If omitted, mpi defaults to "auto".
#     mpi: auto
#     mpi_warning: true
#
# The reserved profile "ase" is read from ASE configuration and is not stored
# here. The reserved "legacy" profile exists only for deprecated options;
# --migrate converts those execution settings to a normal named profile.

# running:
#   empty_spheres: auto
#   print_output: info

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


def _requested_cleanup(args):
    if args.keep_init:
        return "keep"
    if args.clean_init:
        return "clean"
    if args.delete_init:
        return "delete"
    return None


def _migrate(args, legacy_file, file):
    from ...configuration import definition
    from ...configuration_migration import (
        clean_legacy_configuration,
        migrate_legacy_configuration,
    )

    if Path(file).is_file():
        print(f"YAML configuration already exists: '{file}'. Nothing was migrated.")
        return
    if not Path(legacy_file).is_file():
        print(f"No Python configuration hook found at '{legacy_file}'.")
        return

    migrate_legacy_configuration(legacy_file, file, definition)
    print(f"Migrated '{legacy_file}' to '{file}'.")

    action = _requested_cleanup(args) or "clean"
    if action == "keep":
        print(f"Kept Python configuration hook unchanged: '{legacy_file}'.")
    elif action == "clean":
        count = clean_legacy_configuration(legacy_file)
        print(f"Removed {count} migrated configuration assignment(s) from '{legacy_file}'.")
    else:
        Path(legacy_file).unlink()
        print(f"Deleted Python configuration hook: '{legacy_file}'.")


def _parse_set_value(text):
    """Parse a command-line value using YAML scalar/flow syntax."""
    from ruamel.yaml import YAML

    return YAML(typ="safe").load(text)


def _set_configuration_value(config, name, value):
    """Set one CLI value, creating a named profile when necessary."""
    path = name.split(".")
    if path[0] != "profiles":
        config.find(name).set(value)
        config.save_user_preferences()
        return

    if len(path) not in (2, 3) or not path[1]:
        raise KeyError(name)

    profile_name = path[1]
    if profile_name in config.profiles.reserved_names:
        raise ValueError(f"{profile_name!r} is a reserved profile name")

    if len(path) == 2:
        if not isinstance(value, Mapping):
            raise TypeError(
                f"Profile {profile_name!r} must be set from a mapping"
            )
        profile_values = value
    else:
        profile_values = {path[2]: value}

    # Repeated-container staging validates the complete proposed profile and is
    # atomic, so a bad field/value does not leave a half-created profile behind.
    try:
        config.profiles.set(
            {profile_name: profile_values},
            unknown="fail",
            merge=True,
        )
    except KeyError as error:
        raise ValueError(f"Invalid profile option: {error}") from error
    config.save_user_preferences()


def _unset_configuration_value(config, name):
    """Unset one CLI value, or remove one complete named profile."""
    path = name.split(".")
    if len(path) == 2 and path[0] == "profiles" and path[1]:
        config.profiles.remove(path[1])
    else:
        config.find(name).clear()
    config.save_user_preferences()


def _print_profiles(config):
    """Print public execution profiles, marking an explicitly selected default."""
    default = config.default_profile()
    for name in config.profiles:
        suffix = " (default)" if name == default else ""
        print(f"{name}{suffix}")


def _print_executables():
    """Print declarative information about known SPR-KKR executables."""
    from ...sprkkr.executables import executable_registry

    definitions = tuple(executable_registry.values())
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

    cleanup = _requested_cleanup(args)
    if cleanup and not args.migrate:
        raise SystemExit("--keep-init, --clean-init and --delete-init can only be used with --migrate")

    if args.migrate:
        _migrate(args, legacy_file, file)
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
        from ...common.warnings import DataValidityError

        try:
            value = _parse_set_value(args.set[1])
            _set_configuration_value(config, args.set[0], value)
            print(f"Configuration '{args.set[0]}' has been set to '{args.set[1]}'.")
        except KeyError:
            print(f"Unknown configuration name '{args.set[0]}'.")
        except (DataValidityError, TypeError, ValueError) as error:
            print(
                f"Invalid configuration value for '{args.set[0]}': "
                f"{error}"
            )
        except pyparsing.ParseBaseException:
            print(
                f"Suspicious configuration value '{args.set[1]}'. If you are sure, "
                "please edit the configuration file manually."
            )
        run = False
    unset = getattr(args, "unset", None)
    if unset:
        from ...common.warnings import DataValidityError

        try:
            _unset_configuration_value(config, unset)
            print(f"Configuration '{unset}' has been unset.")
        except KeyError:
            print(f"Unknown configuration name '{unset}'.")
        except (DataValidityError, TypeError, ValueError) as error:
            print(f"Cannot unset configuration '{unset}': {error}")
        run = False
    if getattr(args, "list_profiles", False):
        _print_profiles(config)
        run = False
    if getattr(args, "executables", False):
        _print_executables()
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
            print(f"No Python configuration hook found at '{legacy_file}'.")
        run = False
    if args.edit_init_script:
        if not os.path.isfile(legacy_file):
            print(f"No Python configuration hook found at '{legacy_file}'.")
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
            print(f"Python configuration hook found at:\n    {legacy_file}")
            print()
            print("Run 'ase2sprkkr config --migrate' to migrate declarative settings to config.yaml.")
        else:
            print("# No custom user configuration found")


if __name__ == "__main__":
    main(globals())
