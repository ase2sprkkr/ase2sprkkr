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

from ...common.tools import main  # NOQA
from ...bindings.nomad.nomad_api import NomadApi  # NOQA


help='Upload given datas to NOMAD'
description="""\n Uploads a given output file to a NOMAD. You either have to specify your Nomad username and type your password on the prompt.
Or you can authenticate once using -a switch, then the token will be stored. """


def parser(parser):

    nomad = argparse.ArgumentParser(add_help=False)
    nomad.add_argument('--nomad-url', '-N', help='Nomad URL.', type=str, default=NomadApi.default_api_url)
    nomad.add_argument('--password', '-P', help='Nomad password. Since there is command-line history in the most of shells, it is not reccomended to use.', type=str)

    nomad_auth = argparse.ArgumentParser(add_help=False, parents = [ nomad ])
    group = nomad_auth.add_mutually_exclusive_group()
    group.add_argument('--token', '-t', help='Nomad token. Can be also specified in config.nomad.token.', type=str, required=False)
    group.add_argument('--user', '-u', help='Nomad user name. ', type=str, required=False)

    subs = parser.add_subparsers(required=True)

    sub = subs.add_parser('authenticate', help='Retrieve and store Nomad token for further passwordless authentication.', parents=[nomad])
    group = sub.add_mutually_exclusive_group(required=True)
    group.add_argument('user', help='Nomad user name. ', type=str, nargs='?')
    group.add_argument('--delete-credentials', '-d', help='Clear authentication data.', action='store_true')
    sub.add_argument('--print-token', '-o', help='Print the resulting token.', action='store_true')
    sub.add_argument('--do-not-store', '-n', help='Only obtain the token, do not store it.', action='store_true')
    sub.add_argument('--expires', '-e', help='Validity of the token (in days), the default is one year.', type=int, default=365)
    sub.set_defaults(func=authenticate, api=True, token=False)

    sub = subs.add_parser('token', help='Print the (currently stored) Nomad authentication token')
    sub.set_defaults(func=token, api=False, token=False)

    sub = subs.add_parser('upload', help='Upload a file to nomad', parents=[nomad_auth])
    sub.add_argument('output_files', help='The output file to be uploaded. If you specify more files, they will be uploaded in single upload.', nargs='+', type=str)
    sub.set_defaults(func=upload, api=True, token=True)

    sub = subs.add_parser('zip', help='Create Nomad archive')
    sub.add_argument('zip', help='Name of the archive')
    sub.add_argument('output_files', help='The output file to be uploaded. If you specify more files, they will be uploaded in single upload.', nargs='+', type=str)
    sub.set_defaults(func=zipp, api=False, token=False)


# api object used in functions
api=None


def authenticate(args):
    from ase2sprkkr.configuration import config
    if args.delete_credentials:
        token=None
        config.nomad.token.set_permanent(token, r"Authentication token to Nomad \(written by 'ase2sprkkr nomad authenticate [^\s]+'\)", True)
        return
    else:
        token=retrieve_token(args, args.expires * 24 * 3600)
    if not args.do_not_store:
        config.nomad.token.set_permanent(token, f"Authentication token to Nomad (written by 'ase2sprkkr nomad authenticate {args.user}')")

    if args.print_token:
        print(token)


def retrieve_token(args, expires=None):
    import getpass
    password = args.password
    if not password:
        password = getpass.getpass(prompt='Nomad password: ')
    token = api.get_authentication_token(args.user, password, expires=expires)
    return token


def get_token(args):
    from ase2sprkkr.configuration import config
    if args.user:
        token=retrieve_token(args)
    else:
        token=config.nomad.token()
    if not token:
        raise ValueError("No Nomad authentization token: "
                         "please either supply -u parameter, or run `nomad authenticate`")
    return token


def token(args):
    from ase2sprkkr.configuration import config
    print(config.nomad.token())


def gather_files(output_files, name=None):
    from ...bindings.nomad.nomad import NomadArchive
    arch = NomadArchive(name)
    for i in output_files:
        arch.add_entry(i)
    arch.finalize()
    return arch


def upload(args):
    arch = gather_files(args.output_files)
    api.upload(arch.file)


def zipp(args):
    gather_files(args.output_files, args.zip)


def run(args):

    global api

    if args.api:
        api = NomadApi(args.nomad_url)
        if args.token:
            api.token = get_token(args)
    args.func(args)
