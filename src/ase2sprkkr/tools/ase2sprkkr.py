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

if not __package__:
    path = str(Path(__file__).resolve().parents[1])
    sys.path.append(path)
    import os
    spec = importlib.util.spec_from_file_location("ase2sprkkr", os.path.join(path, '__init__.py'))
    ase2sprkkr = importlib.util.module_from_spec(spec)
    sys.modules["ase2sprkkr"] = ase2sprkkr
    spec.loader.exec_module(ase2sprkkr)
    __package__ = 'ase2sprkkr.tools'

import ase2sprkkr.tools.commands as commands # NOQA


def run():

  parser = argparse.ArgumentParser(
      description='ASE2SPRKKR tool: tool for visualising SPRKKR result',
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog='You can install autocompleting for bash and zsh by running/adding to the .bashrc: \n'
             'eval "$(register-python-argcomplete ase2sprkkr)"'
  )
  parser.add_argument('--version', '-v', help='Print the version of ASE2SPRKKR.', action='store_true')
  parser.add_argument('--debug', '-G', help='Raise a debugger on an unhandled exception.', action='store_true')
  parser.add_argument('--profile', '-P', help='Run a python profiler on the command.', action='store_true')

  subs = parser.add_subparsers( dest = 'ase2sprkkr_command', description='Run ase2sprkkr <subcommand> -h for futhrer info')

  names = (i for i in pkgutil.iter_modules(commands.__path__))
  im = importlib.import_module
  modules = ( im(commands.__name__ + '.' + i.name, __package__) for i in names )
  modules = { m.__name__.rsplit('.',1)[1]: m for m in modules if hasattr(m, 'parser') }
  unknowns = {}

  for name, m in modules.items():
      sub = subs.add_parser( name,
                             help=m.help,
                             formatter_class=argparse.RawDescriptionHelpFormatter,
                            description = m.help + '\n' + m.description )
      if hasattr(m, 'unknowns'):
          unknowns[name] = m.unknowns
      m.parser( sub )

  argcomplete.autocomplete(parser)

  args, remainder = parser.parse_known_args()
  if remainder:
      where = unknowns.get(args.ase2sprkkr_command, None)
      if where is None:
          parser.parse_args()
      else:
          where = getattr(args, where)
          where += remainder

  help = True
  if args.debug:
      from ase2sprkkr.common.debug import add_debug_hook
      add_debug_hook()
  del args.debug
  del args.profile

  if args.version:
      import ase2sprkkr.version
      print(ase2sprkkr.version.__version__)
      help=False
  if args.ase2sprkkr_command is None:
      if help:
          parser.print_help()
  else:
      action = modules[ args.ase2sprkkr_command ].run
      del args.ase2sprkkr_command
      action(args)


if __name__ == "__main__":
    if '-P' in sys.argv or '--profile' in sys.argv:
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
