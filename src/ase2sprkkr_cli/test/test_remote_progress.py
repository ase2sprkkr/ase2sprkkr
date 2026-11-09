"""User-facing remote options and SCF progress rendering."""

import argparse
import json
from types import SimpleNamespace

import pytest

from ase2sprkkr_cli.commands import remote


def _parser():
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    return parser


def test_follow_defaults_to_current_position():
    args = _parser().parse_args(["follow", "opaque"])
    assert not args.raw
    assert not args.from_start
    assert not args.info
    assert args.regex is None


def test_follow_options():
    args = _parser().parse_args([
        "follow", "opaque", "--info", "--regex", "WARN", "--regex", "ERR", "--from-start",
    ])
    assert args.info
    assert args.regex == ["WARN", "ERR"]
    assert args.from_start


def test_status_includes_scf_iteration(monkeypatch, capsys):
    from kkrunner.remote import service
    monkeypatch.setattr(service, "status_remote", lambda value: {"state": "running"})
    monkeypatch.setattr(remote, "_output_path", lambda args: "scf.out")
    monkeypatch.setattr(remote, "_scf_progress", lambda path: {
        "iteration": 15, "error": 0.001, "b_error": 0.01,
        "energy": {"EF": 0.5}, "converged": False,
    })
    args = _parser().parse_args(["status", "opaque"])
    remote.run(args, None)
    data = json.loads(capsys.readouterr().out)
    assert data["state"] == "running"
    assert data["scf"]["iteration"] == 15
    assert data["scf"]["error"] == 0.001


def test_status_only_does_not_parse(monkeypatch, capsys):
    from kkrunner.remote import service
    monkeypatch.setattr(service, "status_remote", lambda value: {"state": "running"})
    monkeypatch.setattr(remote, "_scf_progress", lambda path: pytest.fail("unexpected parse"))
    args = _parser().parse_args(["status", "opaque", "--status-only"])
    remote.run(args, None)
    assert json.loads(capsys.readouterr().out) == {"state": "running"}


@pytest.mark.parametrize("action", ["submit", "recover", "ack", "query-files"])
def test_internal_commands_not_exposed(action):
    with pytest.raises(SystemExit):
        _parser().parse_args([action])


def test_follow_regex_decodes_partial_lines(monkeypatch, capsys):
    from kkrunner.remote import service
    calls = []

    def follower(remote_id, stdout, stderr, **kw):
        calls.append(kw)
        stdout.write(b"ignore\nWA")
        stdout.write(b"RN here\n")
        stderr.write(b"ERR here\n")
        return {"state": "finished"}

    monkeypatch.setattr(service, "follow_remote", follower)
    monkeypatch.setattr(remote, "_paths", lambda value: {
        "output_file": "/tmp/scf.out", "error_file": "/tmp/scf.err",
    })
    args = _parser().parse_args(["follow", "opaque", "--regex", "WARN", "--from-start"])
    remote.run(args, None)
    captured = capsys.readouterr()
    assert captured.out == "WARN here\n"
    assert captured.err == ""
    assert calls[0]["stdout_offset"] == 0


def test_follow_info_outputs_incremental_cycle(monkeypatch, capsys):
    from kkrunner.remote import service

    def follower(remote_id, stdout, stderr, **kw):
        stdout.write(b"  8 ERR 0.123E-04 0.456E-06 EF 0.6 0.01 D 12 M 1.2 0.3\n")
        stdout.write(b" ETOT -42.5 SCF - cycle not converged - continue\n")
        return {"state": "finished"}

    monkeypatch.setattr(service, "follow_remote", follower)
    monkeypatch.setattr(remote, "_paths", lambda value: {
        "output_file": None, "error_file": None,
    })
    args = _parser().parse_args(["follow", "opaque", "--info"])
    remote.run(args, None)
    captured = capsys.readouterr()
    assert "Iteration     8 error " in captured.out
    assert "ETOT" not in captured.out
