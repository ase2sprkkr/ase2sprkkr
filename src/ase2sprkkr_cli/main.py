#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
"""Lightweight CLI dispatcher, independent of ase2sprkkr imports."""

import argparse
import importlib
import os
import pkgutil


class ArgparseError(ValueError):
    def __init__(self, parser, message):
        super().__init__(message)
        self.parser = parser


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise ArgparseError(self, message)

    def exit(self, status=0, message=None):
        if status == 0:
            raise SystemExit(0)
        raise ArgparseError(self, message or "")

    def format_help(self):
        result = super().format_help()
        extra_help = getattr(self, "extra_help", None)
        if extra_help is not None:
            extra = extra_help()
            if extra:
                result += "\n" + str(extra).rstrip() + "\n"
        return result


def _plugin_warning(message):
    import sys

    print(f"ase2sprkkr: CLI plugin {message}", file=sys.stderr)


def _plugin_packages():
    """Load packages exposed by installed 'ase2sprkkr:cli' entry points."""
    from importlib.metadata import entry_points

    try:
        entries = entry_points()
        if hasattr(entries, "select"):
            selected = entries.select(group="ase2sprkkr", name="cli")
        else:  # Python 3.9's mapping-style entry_points() result
            selected = (entry for entry in entries.get("ase2sprkkr", ()) if entry.name == "cli")
    except Exception as exc:
        _plugin_warning(f"discovery failed: {exc}")
        return

    # A deterministic order makes conflicts between independently installed
    # plugins predictable; native commands are always registered first.
    for entry in sorted(selected, key=lambda entry: entry.value):
        try:
            package = entry.load()
            if not hasattr(package, "__path__") or not hasattr(package, "__name__"):
                raise TypeError("entry point must expose a Python package")
        except Exception as exc:
            _plugin_warning(f"{entry.value!r} could not be loaded: {exc}")
            continue
        yield package


def run(argv=None):
    parser = Parser(
        description="ASE2SPRKKR tool — command-line interface for ASE2SPRKKR functionalities",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Enable shell completion with: eval "$(register-python-argcomplete ase2sprkkr)"',
    )
    parser.add_argument("--version", "-v", action="store_true", help="Show ASE2SPRKKR version.")
    parser.add_argument("--debug", "-G", action="store_true", help="Enable debugger on unhandled exceptions.")
    parser.add_argument("--profile", "-P", action="store_true", help="Run a Python profiler on the command.")
    subs = parser.add_subparsers(dest="ase2sprkkr_command", description="Use <command> -h for details.")

    from . import commands
    if argv is None:
        import sys
        requested_args = sys.argv[1:]
    else:
        requested_args = argv
    from .errors import CLIUsageError

    modules = {}
    aliases = {}
    parsers = {}
    unknowns = {}
    claimed_names = set()

    def register_command(name, command, *, plugin=False, source=None):
        if not hasattr(command, "parser"):
            return
        try:
            command_aliases = getattr(command, "aliases", None) or []
            if not isinstance(command_aliases, (list, tuple)):
                raise TypeError("aliases must be a list or tuple")
            names = (name, *command_aliases)
            if any(not isinstance(value, str) or not value or value.startswith("-") for value in names):
                raise ValueError("command names and aliases must be nonempty strings without a leading '-'")
            if len(set(names)) != len(names):
                raise ValueError("command name and aliases must be distinct")
            conflicts = claimed_names.intersection(names)
            if conflicts:
                raise ValueError(f"command name/alias conflicts with an existing command: {', '.join(sorted(conflicts))}")
            help_text, description = command.help, command.description
            if not isinstance(help_text, str) or not isinstance(description, str):
                raise TypeError("help and description must be strings")
            if not callable(command.parser) or not callable(command.run):
                raise TypeError("parser and run must be callable")
        except (AttributeError, TypeError, ValueError) as exc:
            if not plugin:
                raise
            _plugin_warning(f"{source}.{name} skipped: {exc}")
            return

        sub = subs.add_parser(
            name, help=help_text, aliases=command_aliases,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=help_text + "\n" + description,
        )
        sub.extra_help = getattr(command, "extra_help", None)
        # Native commands retain their existing lazy_parser setting. Third-party
        # parsers are configured only when requested, keeping unrelated plugin
        # dependencies out of ordinary CLI invocations and root --help.
        if (not plugin and not getattr(command, "lazy_parser", False)) or any(
            value in requested_args for value in names
        ):
            command.parser(sub)
        modules[name] = command
        parsers[name] = sub
        claimed_names.update(names)
        for alias in command_aliases:
            aliases[alias] = name
        if hasattr(command, "unknowns"):
            unknowns[name] = command.unknowns

    for discovered in pkgutil.iter_modules(commands.__path__):
        command = importlib.import_module(f"{commands.__name__}.{discovered.name}")
        register_command(discovered.name.replace("_", "-"), command)

    for package in _plugin_packages():
        try:
            discovered_modules = list(pkgutil.iter_modules(package.__path__))
        except Exception as exc:
            _plugin_warning(f"{package.__name__!r} could not be inspected: {exc}")
            continue
        for discovered in discovered_modules:
            module_name = f"{package.__name__}.{discovered.name}"
            try:
                command = importlib.import_module(module_name)
            except Exception as exc:
                _plugin_warning(f"{module_name!r} could not be imported: {exc}")
                continue
            register_command(
                discovered.name.replace("_", "-"), command,
                plugin=True, source=package.__name__,
            )

    if os.environ.get("_ARGCOMPLETE"):
        import argcomplete

        argcomplete.autocomplete(parser)

    def parser_error(failed, message):
        from colorama import Fore, Style, init

        init()
        print("\n" + Fore.RED + message + Style.RESET_ALL + "\n")
        print("-" * 30 + "\n")
        failed.print_help()
        raise SystemExit(-2)

    try:
        args, remainder = parser.parse_known_args(argv)
    except ArgparseError as exc:
        parser_error(exc.parser, f" ! {exc} !")

    def execute():
        selected = args.ase2sprkkr_command
        selected = aliases.get(selected, selected)
        if remainder:
            target = unknowns.get(selected)
            if target is None:
                parser_error(
                    parsers.get(selected, parser),
                    f" ! ERROR: Unknown argument(s) {', '.join(remainder)!r} for {selected or 'ase2sprkkr'} !",
                )
            getattr(args, target).extend(remainder)

        if args.debug:
            from ase2sprkkr.common.debug import add_debug_hook

            add_debug_hook()

        global_args = {"debug": args.debug, "profile": args.profile}
        if args.version:
            from importlib.metadata import PackageNotFoundError, version

            try:
                print(version("ase2sprkkr"))
            except PackageNotFoundError:
                parser_error(parser, " ! ASE2SPRKKR is not installed ! ")
        if selected is None:
            if not args.version:
                parser.print_help()
            return

        for name in ("ase2sprkkr_command", "version", "debug", "profile"):
            delattr(args, name)
        try:
            modules[selected].run(args, global_args)
        except CLIUsageError as exc:
            parser_error(parsers[selected], f" ! {exc} !")

    if args.profile:
        import cProfile
        import io
        import pstats

        profiler = cProfile.Profile()
        profiler.enable()
        try:
            execute()
        finally:
            profiler.disable()
            stream = io.StringIO()
            pstats.Stats(profiler, stream=stream).sort_stats(pstats.SortKey.CUMULATIVE).print_stats(0.1)
            print(stream.getvalue())
    else:
        execute()


if __name__ == "__main__":
    run()
