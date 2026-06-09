#!/usr/bin/env python
# PYTHON_ARGCOMPLETE_OK
"""
The main ase2sprkkr scripts. See the commands subdir for the available commands.
"""

import argparse
import argcomplete
import sys
import pkgutil
import importlib
from pathlib import Path
import os


def fix_package():
    global __package__

    if not __package__:
        path = str(Path(__file__).resolve().parents[1])
        sys.path.append(path)
        spec = importlib.util.spec_from_file_location("ase2sprkkr", os.path.join(path, "__init__.py"))
        ase2sprkkr = importlib.util.module_from_spec(spec)
        sys.modules["ase2sprkkr"] = ase2sprkkr
        spec.loader.exec_module(ase2sprkkr)
        __package__ = "ase2sprkkr.tools"


def run():


    class ArgparseError(ValueError):
        def __init__(self, parser, msg):
            super().__init__(msg)
            self.parser = parser

    class Parser(argparse.ArgumentParser):
        def error(self, message):
            raise ArgparseError(self, message)

        def exit(self, status=0, message=None):
            raise ArgparseError(self, message or "")

    parser = Parser(
        description="ASE2SPRKKR tool — command-line interface for ASE2SPRKKR functionalities",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="You can install autocompleting for bash and zsh by running/adding to the .bashrc: \n"
        'eval "$(register-python-argcomplete ase2sprkkr)"',
    )
    parser.add_argument("--version", "-v", help="Print the version of ASE2SPRKKR.", action="store_true")
    parser.add_argument("--debug", "-G", help="Raise a debugger on an unhandled exception.", action="store_true")
    parser.add_argument("--profile", "-P", help="Run a python profiler on the command.", action="store_true")
    # parser.add_argument('--no-user-profile', '-U', help='Do not load the user profile file.', action='store_true')

    subs = parser.add_subparsers(
        dest="ase2sprkkr_command", description="Run ase2sprkkr <subcommand> -h for futhrer info",
    )

    # os.environ['ASE2SPRKKR_NO_USER_PROFILE'] = '1'
    fix_package()
    import ase2sprkkr.tools.commands as commands  # NOQA

    names = (i for i in pkgutil.iter_modules(commands.__path__))
    im = importlib.import_module
    modules = (im(commands.__name__ + "." + i.name, __package__) for i in names)
    modules = {m.__name__.rsplit(".", 1)[1]: m for m in modules if hasattr(m, "parser")}
    unknowns = {}

    aliases = {}
    parsers = {}
    for name, m in modules.items():
        name = name.replace("_", "-")
        alias = getattr(m, "aliases", None)
        sub = subs.add_parser(
            name,
            help=m.help,
            aliases=alias or [],
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=m.help + "\n" + m.description,
        )
        if alias:
            for a in alias:
                aliases[a] = name
        if hasattr(m, "unknowns"):
            unknowns[name] = m.unknowns
        m.parser(sub)
        parsers[name] = sub

    argcomplete.autocomplete(parser)

    def parser_error(parser, msg):
        from colorama import init, Fore, Style
        init()
        print()
        print(Fore.RED + msg + Style.RESET_ALL)
        print()
        print('-'*30)
        print()
        parser.print_help()
        print()
        exit(-2)

    try:
        args, remainder = parser.parse_known_args()
    except ArgparseError as e:
        import re
        msg = str(e)

        def enhance():
            match = re.match(r"^argument ase2sprkkr_command: invalid choice: (.*)", msg)
            if match:
                return "ase2sprkkr: unknown command:" + match.group(1)
            if msg.startswith('the following arguments are required:'):
                return f"{e.parser.prog}: " + msg
            return msg

        msg = enhance()
        msg = f" ! {msg} !"
        parser_error(e.parser, msg)

    module = args.ase2sprkkr_command
    if module:
        module = aliases.get(module, module).replace("-", "_")

    if remainder:
        where = unknowns.get(args.ase2sprkkr_command, None)
        if where is None:
             eparser = parsers.get(module, parser)
             msg = f" ! ERROR: Unknown argument(s) '{', '.join(remainder)}' for the '{eparser.prog}' command. !"
             parser_error(eparser, msg)
        else:
            where = getattr(args, where)
            where += remainder

    help = True
    if args.debug:
        from ase2sprkkr.common.debug import add_debug_hook

        add_debug_hook()

    global_args = {}
    for i in ["debug", "profile"]:
        global_args[i] = getattr(args, i)
        delattr(args, i)

    if args.version:
        import ase2sprkkr.version

        print(ase2sprkkr.version.__version__)
        help = False
    else:
        del args.version

    if args.ase2sprkkr_command is None:
        if help:
            parser.print_help()
    else:
        del args.ase2sprkkr_command
        action = modules[module].run
        action(args, global_args)


if __name__ == "__main__":
    if "-P" in sys.argv or "--profile" in sys.argv:
        import cProfile
        import io
        import pstats

        pr = cProfile.Profile()
        pr.enable()
        run()
        pr.disable()
        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats(0.1)
        print(s.getvalue())
    else:
        run()
