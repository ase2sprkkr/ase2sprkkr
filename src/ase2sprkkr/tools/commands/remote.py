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

    # Used both as a diagnostic command and by remote clients that want to
    # report an unavailable ``mpi='auto'`` launcher before submission.
    actions.add_parser("resolve-mpi")


def _print_json(value):
    print(json.dumps(value), flush=True)


def _finalize_failed_worker_start(directory, request, exc):
    """Turn a failed detached-worker start into an explicit terminal record.

    ``submit_remote`` persists QUEUED metadata before starting the detached
    worker. If process creation itself fails, leaving that record queued would
    make recovery believe a computation exists and ``follow`` could wait
    forever. Only rewrite the exact pre-start state; failures after a worker has
    published a pid/state are left to the normal lifecycle machinery.
    """
    from ...execution.remote import _timestamp, _write_state, status_remote

    run_id = request.get("id")
    if not run_id:
        return
    try:
        state = status_remote(directory, run_id)
    except (FileNotFoundError, KeyError, ValueError, json.JSONDecodeError):
        return
    if state.get("state") != "queued" or state.get("pid") is not None:
        return
    _write_state(
        directory,
        run_id,
        state="failed",
        error=repr(exc),
        return_code=None,
        finished_at=_timestamp(),
    )


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
        try:
            result = submit_remote(request)
        except Exception as exc:
            _finalize_failed_worker_start(args.directory, request, exc)
            raise
        _print_json(result)
    elif args.remote_action == "status":
        # A status query is also the recovery point for a target-side PSI/J
        # agent killed by logout/reboot. The persisted native scheduler ids are
        # sufficient for the replacement agent to attach instead of resubmit.
        from ...execution.remote_psij import ensure_remote_psij_agent

        ensure_remote_psij_agent(args.directory, args.run_id)
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
        from ...execution.base import LaunchSpec
        from ...execution.mpi import resolve_mpi_command

        request = json.load(sys.stdin)
        launch = LaunchSpec.from_mpi(request.get("mpi"))
        _print_json({"mpi": resolve_mpi_command(launch.mpi)})
    else:
        raise ValueError("Unknown remote action: {}".format(args.remote_action))
