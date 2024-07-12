#!/usr/bin/env python
"""
The main ase2sprkkr scripts. See the commands subdir for the available commands.
"""
import argparse
import sys
import pkgutil
import importlib
from pathlib import Path

if not __package__:
  __package__ = 'ase2sprkkr.tools'
sys.path.append(str(Path(__file__).resolve().parents[2]))

import ase2sprkkr.tools.commands as commands # NOQA


def run():

  parser = argparse.ArgumentParser(
      description='ASE2SPRKKR tool: tool for visualising SPRKKR result',
      formatter_class=argparse.RawDescriptionHelpFormatter
  )
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

  args, remainder = parser.parse_known_args()
  if remainder:
      where = unknowns.get(args.ase2sprkkr_command, None)
      if where is None:
          parser.parse_args()
      else:
          where = getattr(args, where)
          where += remainder

  if args.ase2sprkkr_command is None:
      parser.print_help()
  else:
      action = modules[ args.ase2sprkkr_command ].run
      del args.ase2sprkkr_command
      action(args)


if __name__ == "__main__":
    run()
