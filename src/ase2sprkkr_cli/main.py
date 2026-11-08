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
    for discovered in pkgutil.iter_modules(commands.__path__):
        command = importlib.import_module(f"{commands.__name__}.{discovered.name}")
        if not hasattr(command, "parser"):
            continue
        name = discovered.name.replace("_", "-")
        modules[name] = command
        command_aliases = getattr(command, "aliases", None) or []
        sub = subs.add_parser(
            name, help=command.help, aliases=command_aliases,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=command.help + "\n" + command.description,
        )
        sub.extra_help = getattr(command, "extra_help", None)
        if not getattr(command, "lazy_parser", False) or name in requested_args:
            command.parser(sub)
        parsers[name] = sub
        for alias in command_aliases:
            aliases[alias] = name
        if hasattr(command, "unknowns"):
            unknowns[name] = command.unknowns

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
