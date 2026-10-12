"""Internal command used by :class:`ase2sprkkr.execution.RemoteExecutor`."""

import argparse
import json
import sys


help = "Control ASE2SPRKKR execution on this machine."
description = "Internal helper for remote ASE2SPRKKR execution."


def parser(parser):
    actions = parser.add_subparsers(dest="remote_action", required=True)

    submit = actions.add_parser("submit")
    submit.add_argument("directory")

    status = actions.add_parser("status")
    status.add_argument("directory")

    cancel = actions.add_parser("cancel")
    cancel.add_argument("directory")

    follow = actions.add_parser("follow")
    follow.add_argument("directory")
    follow.add_argument("--stdout-offset", type=int, default=0)
    follow.add_argument("--stderr-offset", type=int, default=0)

    actions.add_parser("resolve-mpi")


def _print_json(value):
    print(json.dumps(value), flush=True)


def run(args, global_args):
    from ...execution.remote import (
        cancel_remote,
        follow_remote,
        resolve_remote_mpi,
        status_remote,
        submit_remote,
    )

    if args.remote_action == "submit":
        request = json.load(sys.stdin)
        request["directory"] = args.directory
        _print_json(submit_remote(request))
    elif args.remote_action == "status":
        _print_json(status_remote(args.directory))
    elif args.remote_action == "cancel":
        _print_json(cancel_remote(args.directory))
    elif args.remote_action == "follow":
        follow_remote(
            args.directory,
            sys.stdout.buffer,
            sys.stderr.buffer,
            stdout_offset=args.stdout_offset,
            stderr_offset=args.stderr_offset,
        )
    elif args.remote_action == "resolve-mpi":
        request = json.load(sys.stdin)
        _print_json({"mpi": resolve_remote_mpi(request.get("mpi"))})
    else:
        raise ValueError("Unknown remote action: {}".format(args.remote_action))
