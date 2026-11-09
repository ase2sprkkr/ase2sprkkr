"""Human-facing run observation, deliberately independent of KKRunner's CLI.

KKRunner owns execution identity, persistence, and byte transport. This CLI
owns argparse, run selection, SCF interpretation, and text filtering.
"""

import argparse
import json
import os
import re
import sys
import tempfile
from pathlib import Path

help = "Inspect, follow, or cancel persisted runs or explicit backend handles."
description = "Show run status and live output, with SPR-KKR-specific progress reporting."
lazy_parser = True


class _FollowSourceAction(argparse.Action):
    """Reject persisted IDs with direct stream paths in either argument order."""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        if namespace.run_id is not None and (
            namespace.stdout_path is not None or namespace.stderr_path is not None
        ):
            raise argparse.ArgumentError(
                self, "--stdout-path and --stderr-path require --remote-id instead of run_id"
            )


def parser(argument_parser):
    actions = argument_parser.add_subparsers(dest="remote_action", required=True)

    def identity(action, *, run_id_action="store"):
        source = action.add_mutually_exclusive_group(required=True)
        source.add_argument(
            "run_id", nargs="?", action=run_id_action, help="persistent run ID (default)"
        )
        source.add_argument(
            "--remote-id", metavar="HANDLE",
            help="use a low-level backend handle instead (no persisted file metadata)",
        )

    status = actions.add_parser("status", help="show run state and latest SCF progress")
    identity(status)
    status.add_argument(
        "--status-only", action="store_true",
        help="omit SCF progress and show backend state only",
    )

    follow = actions.add_parser("follow", help="observe an existing run's output")
    identity(follow, run_id_action=_FollowSourceAction)
    follow.add_argument(
        "--raw", action="store_true",
        help="send original binary stream frames rather than decoded text",
    )
    follow.add_argument("--info", action="store_true", help="show SCF iteration summaries")
    follow.add_argument(
        "--regex", action="append", metavar="PATTERN",
        help="print matching lines; may be repeated and combined with --info",
    )
    follow.add_argument(
        "--from-start", action="store_true",
        help="replay output from the beginning (default: observe new data only)",
    )
    paths = follow.add_argument_group(
        "Direct backend stream paths",
        "Only with --remote-id; persisted runs use their recorded file locations. Both paths may be supplied.",
    )
    paths.add_argument(
        "--stdout-path", action=_FollowSourceAction, help="backend standard-output file"
    )
    paths.add_argument(
        "--stderr-path", action=_FollowSourceAction, help="backend standard-error file"
    )

    cancel = actions.add_parser("cancel", help="cancel an execution")
    identity(cancel)


def _identity(args):
    return args.remote_id if args.remote_id else args.run_id


def _manager_and_record(run_id):
    from ase2sprkkr.sprkkr.projects import project_execution_manager

    manager = project_execution_manager()
    return manager, manager.store.get(run_id)


def _task_is_scf(record):
    task = record.metadata.get("ase2sprkkr")
    return isinstance(task, dict) and str(task.get("task", "")).upper() == "SCF"


def _scf_summary(progress):
    return {
        "iteration": progress["iteration"],
        "error": progress["error"],
        "b_error": progress["b_error"],
        "fermi_energy": progress["energy"]["EF"],
        "spin_moment": progress["moment"]["spin"],
        "orbital_moment": progress["moment"]["orbital"],
        "converged": progress["converged"],
    }


def _snapshot_progress(manager, record):
    """Read SCF progress from this run's declared stdout, including remote jobs.

    A snapshot is observational and does not update the canonical local output
    or the run record. Copying to a temporary destination bounds local memory.
    """
    if not _task_is_scf(record):
        return None
    output = record.files.output_file
    if output is None:
        return None

    from ase2sprkkr.outputs.readers.scf import ScfProgressParser

    if record.state.is_terminal:
        # Never interpret an arbitrary replacement of a completed run's
        # historical output as that run's SCF progress.
        if not manager.output_file_is_current(record):
            return None
        if os.path.isfile(output):
            return ScfProgressParser.latest_from_file(output)
        return None

    with tempfile.TemporaryDirectory(prefix="ase2sprkkr-scf-progress-") as directory:
        filename = Path(directory) / "progress.out"
        try:
            manager.snapshot_file(record.id, output, destination=filename)
        except FileNotFoundError:
            return None
        return ScfProgressParser.latest_from_file(filename)


def _persisted_status(manager, record):
    """Query backend facts without mutating the client-owned lifecycle record."""
    facts = {
        "run_id": record.id,
        "state": record.state.value,
        "return_code": record.return_code,
        "error": record.message,
        "failure_phase": (
            None if record.failure_phase is None else record.failure_phase.value
        ),
    }
    if record.state.is_terminal or record.execution.backend_id is None:
        return facts

    status, = manager.executor_for(record).status(record)
    facts.update(
        state=status.state.value,
        return_code=status.return_code,
        error=status.message,
        failure_phase=(
            None if status.failure_phase is None else status.failure_phase.value
        ),
    )
    return facts


def _status(args):
    identifier = _identity(args)
    if args.remote_id:
        from kkrunner.remote.service import status_remote
        facts = status_remote(identifier)
    else:
        manager, record = _manager_and_record(identifier)
        facts = _persisted_status(manager, record)
        if not args.status_only:
            progress = _snapshot_progress(manager, record)
            if progress is not None:
                facts["scf"] = _scf_summary(progress)
    print(json.dumps(facts), flush=True)


class _TextSink:
    """Decode independent UTF-8 chunks without confusing flush with EOF."""

    def __init__(self, callback, kind):
        import codecs

        self._decoder = codecs.getincrementaldecoder("utf8")(errors="replace")
        self._callback = callback
        self._kind = kind

    def write(self, data):
        text = self._decoder.decode(data)
        if text:
            self._callback(text, self._kind)

    def flush(self):
        # Transport flush is called on every chunk; not an end-of-stream.
        pass

    def finish(self):
        text = self._decoder.decode(b"", final=True)
        if text:
            self._callback(text, self._kind)


class _FramedSink:
    """Preserve the runner's binary stream protocol without parsing it."""

    def __init__(self, kind, destination):
        self.kind = kind
        self.destination = destination

    def write(self, data):
        from kkrunner.remote.stream import encode_stream_frame

        self.destination.write(encode_stream_frame(self.kind, data))

    def flush(self):
        self.destination.flush()


def _observe(args, observer):
    if args.raw and (args.info or args.regex):
        raise ValueError("--raw cannot be combined with --info or --regex")

    if args.raw:
        from kkrunner.remote.stream import STREAM_MAGIC, encode_terminal_frame

        destination = sys.stdout.buffer
        destination.write(STREAM_MAGIC)
        destination.flush()
        stdout, stderr = _FramedSink("out", destination), _FramedSink("err", destination)
        outcome = observer(stdout, stderr)
        destination.write(encode_terminal_frame(outcome))
        destination.flush()
        return outcome

    from ase2sprkkr.common.output_printing import output_listeners

    listeners, _ = output_listeners(tuple(args.regex or ()) or (False if args.info else True))
    progress = None
    if args.info:
        from ase2sprkkr.outputs.readers.scf import ScfProgressParser

        progress = ScfProgressParser(on_info=lambda message: sys.stdout.write(message))

    def emit(text, kind):
        if kind == "out" and progress is not None:
            progress.feed(text)
        for listener in listeners:
            listener(text, kind)

    stdout, stderr = _TextSink(emit, "out"), _TextSink(emit, "err")
    try:
        return observer(stdout, stderr)
    finally:
        stdout.finish()
        stderr.finish()
        if progress is not None:
            progress.finish()
        for listener in listeners:
            listener.flush()


def _follow(args):
    identifier = _identity(args)
    if args.remote_id:
        from kkrunner.remote.service import follow_remote

        def observe(stdout, stderr):
            # In direct mode no manifest is available; caller may explicitly
            # provide paths for scheduler output files.
            offset = 0 if args.from_start else None
            return follow_remote(
                identifier, stdout, stderr,
                stdout_offset=offset, stderr_offset=offset,
                stdout_path=args.stdout_path, stderr_path=args.stderr_path,
            )
    else:
        manager, _ = _manager_and_record(identifier)

        def observe(stdout, stderr):
            return manager.follow_output(
                identifier, stdout, stderr, from_start=args.from_start,
            )

    return _observe(args, observe)


def _cancel(args):
    identifier = _identity(args)
    if args.remote_id:
        from kkrunner.remote.service import cancel_remote

        facts = cancel_remote(identifier)
    else:
        manager, _ = _manager_and_record(identifier)
        facts = manager.cancel_run(identifier)
    print(json.dumps(facts), flush=True)


def run(args, global_args):
    if args.remote_action == "status":
        return _status(args)
    if args.remote_action == "follow":
        return _follow(args)
    if args.remote_action == "cancel":
        return _cancel(args)
    raise ValueError(f"Unknown remote action: {args.remote_action}")
