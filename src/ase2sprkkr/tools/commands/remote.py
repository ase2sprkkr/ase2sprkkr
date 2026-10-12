"""Control remote ASE2SPRKKR execution and workspace metadata."""

import json
import re
import sys


help = "Control ASE2SPRKKR execution on this machine."
description = "Submit, inspect, follow, cancel and clean remote ASE2SPRKKR runs."


def _duration(value):
    match = re.fullmatch(r"(\d+(?:\.\d+)?)([smhdw])", value.strip().lower())
    if not match:
        raise ValueError("Duration must use s, m, h, d or w, for example 30d")
    amount = float(match.group(1))
    factor = {
        "s": 1,
        "m": 60,
        "h": 3600,
        "d": 86400,
        "w": 604800,
    }[match.group(2)]
    return amount * factor


def parser(parser):
    actions = parser.add_subparsers(dest="remote_action", required=True)

    submit = actions.add_parser("submit")
    submit.add_argument("directory")

    status = actions.add_parser("status")
    status.add_argument("directory")
    status.add_argument("run_id", nargs="+")

    cancel = actions.add_parser("cancel")
    cancel.add_argument("directory")
    cancel.add_argument("run_id")

    follow = actions.add_parser("follow")
    follow.add_argument("directory")
    follow.add_argument("run_id")
    follow.add_argument("--stdout-offset", type=int, default=0)
    follow.add_argument("--stderr-offset", type=int, default=0)

    runs = actions.add_parser("runs", help="List runs recorded in a workspace")
    runs.add_argument("directory")

    clean = actions.add_parser(
        "clean", help="Remove metadata for completed runs in a workspace"
    )
    clean.add_argument("directory")
    clean.add_argument(
        "--older-than",
        type=_duration,
        metavar="AGE",
        help="Only remove completed runs older than AGE, e.g. 30d or 12h",
    )

    # Kept as a diagnostic/compatibility command. Normal remote execution no
    # longer calls it: launch intent is resolved by the target-side executor.
    actions.add_parser("resolve-mpi")


def _print_json(value):
    print(json.dumps(value), flush=True)


def run(args, global_args):
    from ...execution.remote import (
        cancel_remote,
        clean_remote_runs,
        follow_remote,
        list_remote_runs,
        status_remote,
        submit_remote,
    )

    if args.remote_action == "submit":
        request = json.load(sys.stdin)
        request["directory"] = args.directory
        _print_json(submit_remote(request))
    elif args.remote_action == "status":
        states = [status_remote(args.directory, run_id) for run_id in args.run_id]
        _print_json(states[0] if len(states) == 1 else states)
    elif args.remote_action == "cancel":
        _print_json(cancel_remote(args.directory, args.run_id))
    elif args.remote_action == "follow":
        follow_remote(
            args.directory,
            args.run_id,
            sys.stdout.buffer,
            sys.stderr.buffer,
            stdout_offset=args.stdout_offset,
            stderr_offset=args.stderr_offset,
        )
    elif args.remote_action == "runs":
        _print_json(list_remote_runs(args.directory))
    elif args.remote_action == "clean":
        _print_json(clean_remote_runs(args.directory, older_than=args.older_than))
    elif args.remote_action == "resolve-mpi":
        from ...configuration import mpi_runner

        request = json.load(sys.stdin)
        _print_json({"mpi": mpi_runner(request.get("mpi"))})
    else:
        raise ValueError("Unknown remote action: {}".format(args.remote_action))
