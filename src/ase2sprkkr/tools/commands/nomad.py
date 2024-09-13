#!/usr/bin/env python
"""
This is a sctipt to visualise in_struct.inp files. Run it to see the doc.
"""
from pathlib import Path
import sys
import argparse

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ... import InputParameters, OutputFile  # NOQA
from ...common.tools import main  # NOQA


help='Upload given datas to NOMAD'
description="""\n Uploads a given output file to a NOMAD """


def parser(parser):
    parser.add_argument('output_files', help='The output file to be uploaded. If you specify more files, they will be uploaded in single upload', nargs='*', type=str)
    parser.add_argument('--zip', '-z', type=str, help='Create NOMAD archive zip file', default=None)
    parser.add_argument('--token', '-t', help='Nomad token. Can be also specified in config.nomad.token', type=str)
    action = getattr(argparse, 'BooleanOptionalAction', 'store_true')
    parser.add_argument('--verbose', '-v', help='Give a more verbose info. It is set automatically, if the output is limited.', action=action, default=None)


def run(args):

    from ...bindings.nomad.nomad import NomadArchive
    arch = NomadArchive(args.zip)
    for i in args.output_files:
        arch.add_entry(i)
    arch.finalize()
    if arch.zip:
        exit()
