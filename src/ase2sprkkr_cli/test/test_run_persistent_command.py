"""CLI project tracking uses the existing calculator and KKRunner lifecycle."""

import argparse
from types import SimpleNamespace
from uuid import uuid4

import pytest

from ase2sprkkr_cli.commands import run as command
from ase2sprkkr_cli.errors import CLIUsageError


def _args(*argv):
    parser = argparse.ArgumentParser()
    command.parser(parser)
    return parser.parse_args(argv)


def test_cli_accepts_execution_profile_project_id_and_submit_only():
    project_id = uuid4().hex
    args = _args(
        "sample.pot", "--profile", "cluster", "--project-id", project_id,
        "--submit-only", "--directory", "work",
    )
    assert args.profile == "cluster"
    assert args.project_id == project_id
    assert args.submit_only
    assert args.directory == "work"


def test_cli_profile_is_independent_of_global_python_profiler(monkeypatch):
    from ase2sprkkr_cli.main import run as dispatch

    received = []
    monkeypatch.setattr(command, "run", lambda args, global_args: received.append(
        (args.profile, global_args["profile"])
    ))
    dispatch(["run", "sample.pot", "--profile", "cluster"])
    assert received == [("cluster", False)]


def test_submit_only_needs_project_before_loading_calculation():
    with pytest.raises(CLIUsageError, match="requires --project-id"):
        command.run(_args("missing.pot", "--submit-only"), None)


def test_invalid_project_id_has_helpful_error():
    with pytest.raises(CLIUsageError, match="Invalid project ID"):
        command._project_manager("not-a-uuid")


def test_project_id_resolves_existing_application_project(monkeypatch, tmp_path):
    from kkrunner.projects import ProjectRepository
    from ase2sprkkr.sprkkr import projects

    repository = ProjectRepository(tmp_path / "projects")
    project = repository.create("Example", directory=tmp_path / "calculations")
    monkeypatch.setattr(projects, "default_project_repository", lambda: repository)

    manager, resolved = command._project_manager(project.id)
    assert manager.project_id == project.id
    assert resolved == project
    assert manager.store.directory == tmp_path / "projects_runs"

    with pytest.raises(CLIUsageError, match="does not exist"):
        command._project_manager(uuid4().hex)


def _fake_calculator(monkeypatch, *, backend="managed", state=None):
    import ase2sprkkr
    from kkrunner import RunState

    calls = []
    runner = SimpleNamespace(
        id="run-123",
        client=SimpleNamespace(backend=SimpleNamespace(execution_mode=backend)),
        state=RunState.FINISHED if state is None else state,
        return_code=0,
        failure_message=None,
    )
    runner.submit = lambda: calls.append(("submit",))
    runner.wait = lambda: calls.append(("wait",)) or "result"
    runner.close = lambda: calls.append(("close",))

    class Calculator:
        def __init__(self, *, profile=None):
            calls.append(("profile", profile))

        def create_run(self, **kwargs):
            calls.append(("create", kwargs))
            return runner

        def calculate(self, **kwargs):
            calls.append(("calculate", kwargs))
            return "result"

    monkeypatch.setattr(ase2sprkkr, "SPRKKR", Calculator)
    monkeypatch.setattr(ase2sprkkr.Potential, "from_file", lambda _: "potential")
    return calls


def test_project_submit_only_uses_profile_and_project_directory(monkeypatch, capsys):
    calls = _fake_calculator(monkeypatch)
    project_id = uuid4().hex
    manager = SimpleNamespace(project_id=project_id)
    project = SimpleNamespace(directory="/project/work")
    monkeypatch.setattr(command, "_project_manager", lambda identity: (manager, project))

    command.run(
        _args("sample.pot", "--profile", "cluster", "--project-id", project_id, "--submit-only"),
        None,
    )
    assert calls[0] == ("profile", "cluster")
    assert calls[1][0] == "create"
    assert calls[1][1]["project"] is manager
    assert calls[1][1]["directory"] == "/project/work"
    assert [c[0] for c in calls[2:]] == ["submit"]
    assert "Run ID: run-123" in capsys.readouterr().out


def test_project_waits_and_explicit_directory_overrides_project(monkeypatch, capsys):
    calls = _fake_calculator(monkeypatch)
    monkeypatch.setattr(command, "_project_manager", lambda _: (
        object(), SimpleNamespace(directory="/project/work"),
    ))
    command.run(
        _args("sample.pot", "--project-id", uuid4().hex, "--directory", "/another"),
        None,
    )
    assert calls[1][1]["directory"] == "/another"
    assert [c[0] for c in calls[2:]] == ["submit", "wait"]
    assert "Run ID: run-123" in capsys.readouterr().out


def test_submit_only_rejects_direct_backend_and_closes_unstarted_runner(monkeypatch):
    calls = _fake_calculator(monkeypatch, backend="direct")
    monkeypatch.setattr(command, "_project_manager", lambda _: (
        object(), SimpleNamespace(directory=None),
    ))
    with pytest.raises(CLIUsageError, match="managed backend"):
        command.run(_args("sample.pot", "--project-id", uuid4().hex, "--submit-only"), None)
    assert calls[-1] == ("close",)
    assert all(c[0] != "submit" for c in calls)


def test_normal_run_remains_untracked(monkeypatch):
    calls = _fake_calculator(monkeypatch)
    assert command.run(_args("sample.pot", "--profile", "local"), None) == "result"
    assert calls[0] == ("profile", "local")
    assert calls[1][0] == "calculate"
    assert "project" not in calls[1][1]


def test_persistent_run_reports_nonfinished_state(monkeypatch):
    from ase.calculators.calculator import CalculationFailed
    from kkrunner import RunState

    calls = _fake_calculator(monkeypatch, state=RunState.FAILED)
    monkeypatch.setattr(command, "_project_manager", lambda _: (
        object(), SimpleNamespace(directory=None),
    ))
    with pytest.raises(CalculationFailed, match="calculation failed"):
        command.run(_args("sample.pot", "--project-id", uuid4().hex), None)
    assert [c[0] for c in calls[1:]] == ["create", "submit", "wait"]


def test_global_profiler_has_distinct_option_name(capsys):
    from ase2sprkkr_cli.main import run as dispatch

    with pytest.raises(SystemExit) as exc:
        dispatch(["--help"])
    assert exc.value.code == 0
    assert "--cprofile" in capsys.readouterr().out

    with pytest.raises(SystemExit) as exc:
        dispatch(["run", "--help"])
    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert "--profile NAME" in help_text
    assert "--project-id UUID" in help_text
