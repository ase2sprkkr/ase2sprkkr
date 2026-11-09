"""The end-user CLI delegates to the runner's single remote implementation."""
import argparse
import json

import pytest
from kkrunner_cli import main as cli
from kkrunner.remote import service as remote_service
from ase2sprkkr_cli.commands import remote


@pytest.mark.parametrize('action', ['status', 'cancel'])
def test_user_command_delegates_without_protocol_arguments(monkeypatch, capsys, action):
    calls = []
    monkeypatch.setattr(remote_service, action + '_remote',
                        lambda handle: calls.append(handle) or {'state': 'running'})
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    if action == 'status':
        monkeypatch.setattr(remote, '_output_path', lambda args: None)
    remote.run(parser.parse_args([action, 'opaque-handle']), None)
    assert calls == ['opaque-handle']
    assert json.loads(capsys.readouterr().out) == {'state': 'running'}


def test_follow_exposes_user_stream_options():
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    args = parser.parse_args(['follow', 'opaque', '--live'])
    assert args.live
    assert args.remote_id == 'opaque'


@pytest.mark.parametrize('action', ['submit', 'recover', 'ack', 'ack-submit', 'resolve-mpi'])
def test_internal_protocol_commands_belong_to_runner(action):
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    with pytest.raises(SystemExit):
        parser.parse_args([action])


def test_status_includes_last_completed_scf_iteration(monkeypatch, tmp_path, capsys):
    from ase2sprkkr.outputs.readers.scf import ScfProgressParser
    output = tmp_path / "scf.out"
    output.write_text(
        "TASK = SCF\n"
        "  8 ERR 0.123E-04 0.456E-06 EF 0.6 0.01 D 12 M 1.2 0.3\n"
        " ETOT -42.5 SCF - cycle not converged - continue\n"
        "  9 ERR 0.100E-05 0.000E+00 EF 0.6 0.0 D 12 M 1.2 0.3\n"
    )
    monkeypatch.setattr(remote_service, "status_remote", lambda _: {"state": "running"})
    monkeypatch.setattr(remote_service, "remote_stream_paths",
                        lambda _: {"output_file": str(output), "error_file": None})
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    remote.run(parser.parse_args(["status", "handle"]), None)
    status = json.loads(capsys.readouterr().out)
    assert status["state"] == "running"
    assert status["scf"]["iteration"] == 8
    assert status["scf"]["error"] == pytest.approx(0.123e-4)
    assert status["scf"]["converged"] is False

    remote.run(parser.parse_args(["status", "handle", "--status-only"]), None)
    assert json.loads(capsys.readouterr().out) == {"state": "running"}


def test_follow_live_default_and_from_start(monkeypatch, capsys):
    calls = []
    monkeypatch.setattr(remote_service, "remote_stream_paths",
                        lambda _: {"output_file": None, "error_file": None})
    def following(remote_id, stdout, stderr, **kwargs):
        calls.append(kwargs)
        stdout.write(b"output\n")
        stderr.write(b"warning\n")
        return {"state": "finished"}
    monkeypatch.setattr(remote_service, "follow_remote", following)
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    remote.run(parser.parse_args(["follow", "handle"]), None)
    out = capsys.readouterr()
    assert out.out == "output\n"
    assert out.err == "warning\n"
    assert calls[-1]["stdout_offset"] is None
    remote.run(parser.parse_args(["follow", "handle", "--from-start"]), None)
    assert calls[-1]["stdout_offset"] == 0


def test_follow_info_and_regex_are_additive(monkeypatch, capsys):
    monkeypatch.setattr(remote_service, "remote_stream_paths",
                        lambda _: {"output_file": None, "error_file": None})
    def following(_, stdout, stderr, **kwargs):
        stdout.write(
            b"  8 ERR 0.123E-04 0.456E-06 EF 0.6 0.01 D 12 M 1.2 0.3\n"
            b" ETOT -42.5 SCF - cycle not converged - continue\n"
            b"ERROR: supplemental\n"
        )
        return {"state": "finished"}
    monkeypatch.setattr(remote_service, "follow_remote", following)
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    remote.run(parser.parse_args(["follow", "handle", "--info", "--regex", "ERROR"]), None)
    out = capsys.readouterr().out
    assert "Iteration     8 error" in out
    assert "ERROR: supplemental" in out
    assert "ETOT" not in out


def test_follow_raw_preserves_binary_framing(monkeypatch):
    calls = []
    from kkrunner_cli import main as runner_cli
    monkeypatch.setattr(runner_cli, "run", lambda args, global_args: calls.append(args))
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    remote.run(parser.parse_args(["follow", "handle", "--raw"]), None)
    assert calls[0].framed is True
    assert calls[0].live is True


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
