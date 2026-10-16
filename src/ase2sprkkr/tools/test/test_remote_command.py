import io
from types import SimpleNamespace

from ase2sprkkr.tools.commands import remote as remote_command


def test_remote_status_forwards_only_remote_id(monkeypatch):
    calls = []
    output = []
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_service.status_remote",
        lambda remote_id: calls.append(remote_id) or {"state": "running"},
    )
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(
        SimpleNamespace(remote_action="status", remote_id="opaque-id"),
        None,
    )

    assert calls == ["opaque-id"]
    assert output == [{"state": "running"}]


def test_remote_recovery_reports_unsupported_backend(monkeypatch):
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_service.recover_remote",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            NotImplementedError("no recovery")
        ),
    )
    monkeypatch.setattr(
        remote_command.sys,
        "stdin",
        io.StringIO('{"executor": "local", "submitted_at": "now"}'),
    )
    output = []
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(
        SimpleNamespace(
            remote_action="recover",
            directory="/work",
            run_id="run-1",
        ),
        None,
    )

    assert output == [{"unsupported": True, "error": "no recovery"}]


def test_resolve_mpi_command_preserves_integer_process_count(monkeypatch):
    monkeypatch.setattr(
        "ase2sprkkr.execution.mpi.find_default_mpi_launcher",
        lambda: ("mpirun",),
    )
    monkeypatch.setattr(remote_command.sys, "stdin", io.StringIO('{"mpi": 4}'))
    output = []
    monkeypatch.setattr(remote_command, "_print_json", output.append)

    remote_command.run(SimpleNamespace(remote_action="resolve-mpi"), None)

    assert output == [{"mpi": ("mpirun", "-np", "4")}]
