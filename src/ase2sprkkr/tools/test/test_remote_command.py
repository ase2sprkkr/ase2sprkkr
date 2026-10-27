import argparse
import io
import json
from types import SimpleNamespace

import pytest

import ase2sprkkr.execution.remote_service as remote_service
from ase2sprkkr.execution.remote_protocol import CONTROL_PROTOCOL_VERSION
from ase2sprkkr.tools.commands import remote as remote_command


def remote_args(**kwargs):
    return SimpleNamespace(protocol_version=CONTROL_PROTOCOL_VERSION, **kwargs)


def test_status_command_forwards_opaque_remote_handle(monkeypatch):
    remote_id = "opaque-handle"
    calls = []
    monkeypatch.setattr(
        remote_service,
        "status_remote",
        lambda value: calls.append(value) or {"state": "running"},
    )
    output = []
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(
        remote_args(remote_action="status", remote_id=remote_id),
        None,
    )

    assert calls == [remote_id]
    assert output == [{"state": "running"}]


def test_submit_command_passes_complete_request_to_stateless_service(monkeypatch):
    calls = []
    monkeypatch.setattr(
        remote_service,
        "submit_remote",
        lambda request: calls.append(dict(request))
        or {"remote_id": "opaque", "state": "running"},
    )
    monkeypatch.setattr(
        remote_command.sys,
        "stdin",
        io.StringIO(json.dumps({"id": "run-1", "argv": ["true"]})),
    )
    output = []
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(
        remote_args(remote_action="submit", directory="/remote/work"),
        None,
    )

    assert calls == [
        {"id": "run-1", "argv": ["true"], "directory": "/remote/work"}
    ]
    assert output == [{"remote_id": "opaque", "state": "running"}]


def test_resolve_mpi_command_preserves_integer_process_count(monkeypatch):
    monkeypatch.setattr(
        "ase2sprkkr.execution.mpi.find_default_mpi_launcher",
        lambda: ("mpirun",),
    )
    monkeypatch.setattr(remote_command.sys, "stdin", io.StringIO('{"mpi": 4}'))
    output = []
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(remote_args(remote_action="resolve-mpi"), None)

    assert output == [{"mpi": ("mpirun", "-np", "4")}]


def test_incompatible_protocol_rejected_before_request_is_read(monkeypatch):
    monkeypatch.setattr(remote_command.sys, "stdin", io.StringIO("not JSON"))
    args = SimpleNamespace(remote_action="submit", protocol_version=CONTROL_PROTOCOL_VERSION + 1)
    with pytest.raises(ValueError, match="Unsupported remote control protocol"):
        remote_command.run(args, None)


def test_control_protocol_version_is_required_by_parser():
    parser = argparse.ArgumentParser()
    remote_command.parser(parser)
    with pytest.raises(SystemExit):
        parser.parse_args(["status", "opaque"])
    parsed = parser.parse_args([
        "--protocol-version", str(CONTROL_PROTOCOL_VERSION), "status", "opaque",
    ])
    assert parsed.protocol_version == CONTROL_PROTOCOL_VERSION
