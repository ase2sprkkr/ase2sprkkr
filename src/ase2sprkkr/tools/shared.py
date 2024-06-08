"""
Cli subcommands can be runned on its own.
This method handles it.
"""

import argparse


def main(local):
    parser = argparse.ArgumentParser(
      description=local['description'],
      formatter_class=argparse.RawDescriptionHelpFormatter
    )
    local['parser'](parser)
    args = parser.parse_args()
    local['run'](args)
