"""ASE2SPRKKR owns its CLI grammar and resolves persistent run metadata itself."""

import argparse
import io
import json
from types import SimpleNamespace

import pytest

from ase2sprkkr_cli.commands import remote


def _parse(*args):
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    return parser.parse_args(args)


def test_command_grammar_is_separate_from_kkrunner():
    args = _parse("status", "run-123")
    assert args.run_id == "run-123"
    assert args.remote_id is None
    direct = _parse("status", "--remote-id", "opaque")
    assert direct.run_id is None
    assert direct.remote_id == "opaque"
    assert _parse("follow", "run-123", "--info", "--from-start").info
    assert _parse("follow", "run-123", "--regex", "WARN", "--regex", "ERR").regex == ["WARN", "ERR"]
    assert _parse("cancel", "--remote-id", "opaque").remote_id == "opaque"
    for action in ("submit", "recover", "ack", "ack-submit", "resolve-mpi"):
        with pytest.raises(SystemExit):
            _parse(action)


@pytest.mark.parametrize("action", ["status", "follow", "cancel"])
@pytest.mark.parametrize("identity", [(), ("run-123", "--remote-id", "opaque")])
def test_parser_requires_exactly_one_identity(action, identity):
    with pytest.raises(SystemExit) as error:
        _parse(action, *identity)
    assert error.value.code == 2


@pytest.mark.parametrize("paths", [
    ("--stdout-path", "stdout.out"),
    ("--stderr-path", "stderr.out"),
    ("--stdout-path", "stdout.out", "--stderr-path", "stderr.out"),
])
@pytest.mark.parametrize("paths_first", [False, True])
def test_stream_paths_require_direct_identity_in_either_argument_order(paths, paths_first):
    arguments = (*paths, "run-123") if paths_first else ("run-123", *paths)
    with pytest.raises(SystemExit) as error:
        _parse("follow", *arguments)
    assert error.value.code == 2


@pytest.mark.parametrize("paths_first", [False, True])
def test_direct_follow_accepts_both_stream_paths(monkeypatch, paths_first):
    from kkrunner.remote import service

    calls = []
    monkeypatch.setattr(service, "follow_remote", lambda *args, **kwargs: calls.append(kwargs))
    paths = ("--stdout-path", "stdout.out", "--stderr-path", "stderr.out")
    identity = ("--remote-id", "opaque")
    arguments = (*paths, *identity) if paths_first else (*identity, *paths)
    remote.run(_parse("follow", *arguments), None)
    assert calls[0]["stdout_path"] == "stdout.out"
    assert calls[0]["stderr_path"] == "stderr.out"


def test_invalid_stream_source_is_a_cli_usage_error(capsys):
    from ase2sprkkr_cli.main import run

    with pytest.raises(SystemExit) as error:
        run(["remote", "follow", "run-123", "--stdout-path", "stdout.out"])
    assert error.value.code == -2
    assert "require --remote-id instead of run_id" in capsys.readouterr().out


def _record(*, state="running", filename=None, task="SCF"):
    return SimpleNamespace(
        id="run-123",
        state=SimpleNamespace(value=state, is_terminal=(state == "finished")),
        return_code=None, message=None, failure_phase=None,
        metadata={"ase2sprkkr": {"task": task}},
        files=SimpleNamespace(output_file=filename),
        execution=SimpleNamespace(backend_id="opaque"),
    )


def test_status_uses_persistent_run_and_adds_scf_progress(monkeypatch, capsys):
    record = _record(filename="scf.out")
    status = SimpleNamespace(
        state=SimpleNamespace(value="running"), return_code=None,
        message=None, failure_phase=None,
    )
    manager = SimpleNamespace(
        store=SimpleNamespace(get=lambda run: record),
        executor_for=lambda run: SimpleNamespace(status=lambda *records: (status,)),
    )
    monkeypatch.setattr(remote, "_manager_and_record", lambda run: (manager, record))
    monkeypatch.setattr(remote, "_snapshot_progress", lambda m, r: {
        "iteration": 8, "error": 1e-4, "b_error": 0,
        "energy": {"EF": 0.6},
        "moment": {"spin": 1.2, "orbital": 0.3},
        "converged": False,
    })
    remote.run(_parse("status", "run-123"), None)
    facts = json.loads(capsys.readouterr().out)
    assert facts["run_id"] == "run-123"
    assert facts["state"] == "running"
    assert facts["scf"]["iteration"] == 8


def test_status_only_does_not_snapshot(monkeypatch, capsys):
    record = _record(state="finished")
    monkeypatch.setattr(remote, "_manager_and_record", lambda run: (None, record))
    monkeypatch.setattr(remote, "_snapshot_progress", lambda *args: pytest.fail("unexpected snapshot"))
    remote.run(_parse("status", "run-123", "--status-only"), None)
    assert "scf" not in json.loads(capsys.readouterr().out)


def test_direct_status_avoids_persistence(monkeypatch, capsys):
    from kkrunner.remote import service
    monkeypatch.setattr(remote, "_manager_and_record", lambda _: pytest.fail("unexpected store"))
    monkeypatch.setattr(service, "status_remote", lambda _: {"state": "running"})
    remote.run(_parse("status", "--remote-id", "opaque"), None)
    assert json.loads(capsys.readouterr().out) == {"state": "running"}


def test_follow_uses_manager_observation_and_filters(monkeypatch, capsys):
    seen = []

    def following(run_id, out, err, *, from_start):
        seen.append((run_id, from_start))
        out.write(b"ignore\nWA")
        out.flush()
        out.write(b"RN here\n")
        err.write(b"ERR here\n")
        return {"state": "finished"}

    manager = SimpleNamespace(follow_output=following)
    monkeypatch.setattr(remote, "_manager_and_record", lambda run: (manager, _record()))
    remote.run(_parse("follow", "run-123", "--regex", "WARN", "--from-start"), None)
    captured = capsys.readouterr()
    assert seen == [("run-123", True)]
    assert captured.out == "WARN here\n"
    assert captured.err == ""


def test_follow_info_is_additive_with_regex(monkeypatch, capsys):
    def following(_, out, err, *, from_start):
        out.write(
            b"  8 ERR 0.123E-04 0.456E-06 EF 0.6 0.01 D 12 M 1.2 0.3\n"
            b" ETOT -42.5 SCF - cycle not converged - continue\n"
            b"ERROR: supplemental\n"
        )
        return {"state": "finished"}

    monkeypatch.setattr(remote, "_manager_and_record", lambda run: (
        SimpleNamespace(follow_output=following), _record(),
    ))
    remote.run(_parse("follow", "run-123", "--info", "--regex", "ERROR"), None)
    output = capsys.readouterr().out
    assert "Iteration     8 error" in output
    assert "ERROR: supplemental" in output
    assert "ETOT" not in output


def test_follow_raw_has_framing(monkeypatch):
    from kkrunner.remote.stream import STREAM_MAGIC
    capture = io.BytesIO()
    proxy = SimpleNamespace(buffer=capture, write=lambda *args: None, flush=lambda: None)
    monkeypatch.setattr(remote.sys, "stdout", proxy)

    def following(_, out, err, *, from_start):
        out.write(b"raw")
        err.write(b"error")
        return {"state": "finished"}

    monkeypatch.setattr(remote, "_manager_and_record", lambda run: (
        SimpleNamespace(follow_output=following), _record(),
    ))
    remote.run(_parse("follow", "run-123", "--raw"), None)
    assert capture.getvalue().startswith(STREAM_MAGIC)


def test_cancel_persisted_and_direct(monkeypatch, capsys):
    from kkrunner.remote import service
    monkeypatch.setattr(remote, "_manager_and_record", lambda run: (
        SimpleNamespace(cancel_run=lambda rid: {"state": "running"}), _record(),
    ))
    monkeypatch.setattr(service, "cancel_remote", lambda _: {"state": "cancelled"})
    remote.run(_parse("cancel", "run-123"), None)
    assert json.loads(capsys.readouterr().out)["state"] == "running"
    remote.run(_parse("cancel", "--remote-id", "opaque"), None)
    assert json.loads(capsys.readouterr().out)["state"] == "cancelled"


def test_parser_loads_without_calculation_or_execution_library():
    import subprocess
    import sys

    code = """\
import sys
import argparse
from ase2sprkkr_cli.commands import remote
p = argparse.ArgumentParser()
remote.parser(p)
p.parse_args(["status", "opaque"])
assert "ase2sprkkr" not in sys.modules
assert "kkrunner.remote.service" not in sys.modules
assert "kkrunner.runs.manager" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
