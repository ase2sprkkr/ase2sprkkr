"""Human-facing remote execution commands.

The execution protocol belongs to kkrunner. This module only renders its
output for humans, and deliberately does not import a runner on top-level help.
"""

import json
import os
import re
import sys
from pathlib import Path

help = "Inspect, follow or cancel an execution handle on this host."
description = "Human-readable execution status and output, with raw protocol access."
lazy_parser = True


def _runner_cli():
    try:
        from kkrunner_cli import main
    except ImportError:
        from kkrunner import cli as main
    return main


def parser(argument_parser):
    commands = _runner_cli().user_parser(argument_parser)
    status = commands["status"] if isinstance(commands, dict) else None
    follow = commands["follow"] if isinstance(commands, dict) else None
    if status is None or follow is None:
        # Older kkrunner.cli.user_parser returns None. argparse subparsers
        # expose the actual parsers without rebuilding their definitions.
        group = next(action for action in argument_parser._actions
                     if hasattr(action, "choices") and "status" in action.choices)
        status, follow = group.choices["status"], group.choices["follow"]
    status.add_argument("--status-only", action="store_true",
                        help="Only report backend state, without SCF iteration progress")
    follow.add_argument("--raw", action="store_true",
                        help="Write binary framed stream protocol (including terminal frame)")
    follow.add_argument("--info", action="store_true",
                        help="Print parsed SCF iteration summaries")
    follow.add_argument("--regex", action="append", metavar="PATTERN",
                        help="Print matching output lines; may be repeated")
    follow.add_argument("--from-start", action="store_true",
                        help="Replay saved output from the beginning (otherwise follow live bytes)")
    return commands


def _stream_paths(remote_id, stdout_path=None, stderr_path=None):
    from kkrunner.remote import service
    hint = getattr(service, "remote_stream_paths", None)
    try:
        locations = hint(remote_id) if hint else {}
    except ValueError:
        # An invalid handle will still be rejected by the actual status/follow
        # operation; missing optional metadata must not mask that error.
        locations = {}
    return (stdout_path or locations.get("output_file"),
            stderr_path or locations.get("error_file"))


def _is_scf_output(path):
    if path is None:
        return False
    try:
        with open(path, "rb") as stream:
            header = stream.read(65536)
    except (OSError, ValueError):
        return False
    return bool(re.search(rb"(?m)^\s*TASK\s*=\s*SCF\s*$", header))


def _scf_progress(path):
    if not _is_scf_output(path):
        return None
    from ase2sprkkr.outputs.readers.scf import ScfProgressParser
    try:
        return ScfProgressParser.latest_from_file(path)
    except (OSError, ValueError):
        return None


def _status(args):
    from kkrunner.remote.service import status_remote
    facts = status_remote(args.remote_id)
    if not getattr(args, "status_only", False):
        path, _ = _stream_paths(args.remote_id)
        if path is not None:
            last = _scf_progress(path)
            if last is not None:
                facts = dict(facts, scf={
                    "iteration": last["iteration"],
                    "error": last["error"],
                    "b_error": last["b_error"],
                    "fermi_energy": last["energy"]["EF"],
                    "spin_moment": last["moment"]["spin"],
                    "orbital_moment": last["moment"]["orbital"],
                    "converged": last["converged"],
                })
    print(json.dumps(facts), flush=True)


class _TextSink:
    """Decode stream boundaries safely and apply human-facing output listeners."""

    def __init__(self, listener, kind):
        import codecs
        self._decode = codecs.getincrementaldecoder("utf8")(errors="replace")
        self.listener = listener
        self.kind = kind

    def write(self, data):
        decoded = self._decode.decode(data)
        if decoded:
            self.listener(decoded, self.kind)

    def flush(self):
        # KKRunner invokes flush after each chunk. Do not finalize a decoder
        # while a multibyte character may still be arriving.
        pass

    def finish(self):
        remaining = self._decode.decode(b"", final=True)
        if remaining:
            self.listener(remaining, self.kind)


def _follow(args, global_args):
    if args.raw:
        if args.info or args.regex:
            raise ValueError("--raw cannot be combined with --info or --regex")
        # Delegate binary framing to the protocol owner; no decoding or filtering.
        args.framed = True
        args.live = not args.from_start
        args.stdout_path, args.stderr_path = _stream_paths(
            args.remote_id, args.stdout_path, args.stderr_path
        )
        return _runner_cli().run(args, global_args)

    from ase2sprkkr.common.output_printing import output_listeners
    from kkrunner.remote.service import follow_remote
    out_path, err_path = _stream_paths(args.remote_id, args.stdout_path, args.stderr_path)
    regexes = tuple(args.regex or ())
    if args.info:
        from ase2sprkkr.outputs.readers.scf import ScfProgressParser
        info = ScfProgressParser()
    else:
        info = None

    if args.info or regexes:
        listeners, _ = output_listeners(regexes if regexes else False)
    else:
        listeners, _ = output_listeners(True)

    if info is not None:
        info.on_info = lambda message: sys.stdout.write(message)

    def emit(text, kind):
        if info is not None and kind == "out":
            info.feed(text)
        for listener in listeners:
            listener(text, kind)

    stdout = _TextSink(emit, "out")
    stderr = _TextSink(emit, "err")
    live = not args.from_start
    try:
        return follow_remote(
            args.remote_id, stdout, stderr,
            stdout_offset=None if live or getattr(args, "live_stdout", False)
                          else getattr(args, "stdout_offset", 0),
            stderr_offset=None if live or getattr(args, "live_stderr", False)
                          else getattr(args, "stderr_offset", 0),
            stdout_path=out_path, stderr_path=err_path,
        )
    finally:
        stdout.finish()
        stderr.finish()
        if info is not None:
            info.finish()
        for listener in listeners:
            listener.flush()


def run(args, global_args):
    if args.remote_action == "status":
        return _status(args)
    if args.remote_action == "follow":
        return _follow(args, global_args)
    return _runner_cli().run(args, global_args)
