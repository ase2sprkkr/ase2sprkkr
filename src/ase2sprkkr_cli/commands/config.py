#!/usr/bin/env python
"""Tool for ASE2SPRKKR user configuration."""

from pathlib import Path
from collections.abc import Mapping




help = "Show or edit the current user configuration."
description = "The user configuration is stored in ~/.config/ase2sprkkr/config.yaml on typical Unix/Linux systems."


def parser(parser):
    parser.add_argument("-p", "--path", help="Print the path to config.yaml.", action="store_true")
    parser.add_argument("-P", "--print", help="Print config.yaml.", action="store_true")
    parser.add_argument("-s", "--show", help="Print the resolved configuration.", action="store_true")
    parser.add_argument("-i", "--info", help="Show the description of the configuration options.", action="store_true")
    parser.add_argument(
        "--system",
        action="store_true",
        help=(
            "Permit modifications to ASE configuration files under /etc "
            "when setting or unsetting profiles.ase (or profiles.default "
            "when it resolves to ASE). Otherwise such files are never edited."
        ),
    )
    parser.add_argument(
        "-S",
        "--set",
        nargs=2,
        help=(
            "Set a configuration value in its owning configuration file. "
            "Examples: '-S profiles.default.executable_dir /opt/sprkkr', "
            "'-S profiles.ase.executable_suffix _v8', or "
            "'-S profiles.cluster.mpi 8'."
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
    add_migration_cleanup_arguments(parser)
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


def add_migration_cleanup_arguments(parser):
    """Share migration safety options between config --migrate and migrate."""
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
#     # SSH profiles may choose backend: direct, kkrunner or psij.
#
# "default" is a virtual alias for the selected profile and is not stored.
# The editable "ase" view is backed by ASE's own config.ini, not our YAML.
# The reserved "legacy" profile exists only for deprecated options;
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
    from ase2sprkkr.configuration import definition
    from ase2sprkkr.configuration_migration import (
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


def _set_configuration_value(config, name, value, *, system=False):
    """Persist a configuration value to the correct owning configuration."""
    path = name.split(".")
    if path[0] != "profiles":
        if system:
            raise ValueError("--system is only valid for ASE profile updates")
        config.find(name).set(value)
        config.save_user_preferences()
        return

    if len(path) not in (2, 3) or not path[1]:
        raise KeyError(name)
    profile_name = path[1]
    if profile_name == "default":
        profile_name = config.profiles.default_write_target()

    if len(path) == 2:
        if not isinstance(value, Mapping):
            raise TypeError(f"Profile {profile_name!r} must be set from a mapping")
        profile_values = value
    else:
        profile_values = {path[2]: value}

    if profile_name == "ase":
        # ASE is a virtual profile backed by ASE's INI file, not our YAML.
        profile = config.profiles["ase"]
        profile.set(profile_values)
        profile.store_permanent(system=system)
        return
    if system:
        raise ValueError("--system is only valid for ASE profile updates")
    if profile_name in config.profiles.reserved_names:
        raise ValueError(f"{profile_name!r} is a reserved profile name")

    try:
        config.profiles.set(
            {profile_name: profile_values}, unknown="fail", merge=True,
        )
    except KeyError as error:
        raise ValueError(f"Invalid profile option: {error}") from error
    config.save_user_preferences()


def _unset_configuration_value(config, name, *, system=False):
    """Remove a setting from its owning configuration, without copying it."""
    path = name.split(".")
    if path[0] != "profiles":
        if system:
            raise ValueError("--system is only valid for ASE profile updates")
        config.find(name).clear()
        config.save_user_preferences()
        return

    if len(path) not in (2, 3) or not path[1]:
        raise KeyError(name)
    profile_name = path[1]
    if profile_name == "default":
        if len(path) != 3:
            raise ValueError("'default' is an alias; unset a field, not the profile")
        profile_name = config.profiles.default_write_target()

    if profile_name == "ase":
        if len(path) != 3:
            raise ValueError("'ase' is a virtual profile; unset a field instead")
        profile = config.profiles["ase"]
        profile.unset(path[2])
        profile.store_permanent(system=system)
        return
    if system:
        raise ValueError("--system is only valid for ASE profile updates")
    if len(path) == 2:
        config.profiles.remove(profile_name)
    else:
        config.find(f"profiles.{profile_name}.{path[2]}").clear()
    config.save_user_preferences()


def _print_profiles(config):
    """Print public execution profiles, marking the resolved default."""
    default = config.profiles.resolve_name()
    for name in config.profiles:
        suffix = " (default)" if name == default else ""
        print(f"{name}{suffix}")


def run(args, global_args):
    import os
    import pyparsing
    from ase2sprkkr.configuration import (
        config,
        legacy_user_preferences_file,
        user_preferences_file,
    )

    file = user_preferences_file()
    legacy_file = legacy_user_preferences_file()
    run = True

    cleanup = _requested_cleanup(args)
    if args.system and not (args.set or args.unset):
        raise SystemExit("--system requires --set or --unset on an ASE profile")
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
        from ase2sprkkr.common.warnings import DataValidityError

        try:
            value = _parse_set_value(args.set[1])
            _set_configuration_value(config, args.set[0], value, system=args.system)
            print(f"Configuration '{args.set[0]}' has been set to '{args.set[1]}'.")
        except KeyError:
            print(f"Unknown configuration name '{args.set[0]}'.")
        except (DataValidityError, TypeError, ValueError, PermissionError) as error:
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
        from ase2sprkkr.common.warnings import DataValidityError

        try:
            _unset_configuration_value(config, unset, system=args.system)
            print(f"Configuration '{unset}' has been unset.")
        except KeyError:
            print(f"Unknown configuration name '{unset}'.")
        except (DataValidityError, TypeError, ValueError, PermissionError) as error:
            print(f"Cannot unset configuration '{unset}': {error}")
        run = False
    if getattr(args, "list_profiles", False):
        _print_profiles(config)
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

