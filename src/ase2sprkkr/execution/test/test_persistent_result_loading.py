"""Historical SPR-KKR results are reconstructed by the backend, not the GUI."""

from __future__ import annotations

import pytest

from ase2sprkkr.execution import (
    ExecutionManager, FileEntry, FileManifest, ProjectRepository, RunState, RunStore,
)
from ase2sprkkr.outputs.reader_registry import reader_registry
from ase2sprkkr.outputs.task_result import TaskResult


def _stored_run(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Test")
    manager = ExecutionManager(
        store=RunStore(tmp_path / "runs"), project_repository=projects, project=project
    )
    output = tmp_path / "saved.out"
    output.write_text("test output")
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    manager.store.save(record)
    record = manager.complete(
        record, state=RunState.FINISHED, return_code=0, artifacts_available=True
    )
    return manager, record, output


def test_load_result_uses_persisted_task_and_attaches_record(monkeypatch, tmp_path):
    manager, record, output = _stored_run(tmp_path)
    calls = []

    class Reader:
        def __init__(self, input_parameters, calculator, directory, task_name=None):
            calls.append(("reader", task_name, str(directory)))

        def read_from_file(self, filename, error=None, return_code=0, input_file=None):
            calls.append(("read", str(filename), return_code))
            return type("Result", (), {})()

    monkeypatch.setattr(
        reader_registry, "class_for_task_name",
        lambda name: calls.append(("task", name)) or Reader,
    )

    result = manager.load_result(record.id)
    assert result.run_record.id == record.id
    assert calls[0] == ("task", "DOS")
    assert calls[-1] == ("read", str(output), 0)


def test_load_result_rejects_modified_historical_output(tmp_path):
    manager, record, output = _stored_run(tmp_path)
    output.write_text("overwritten")
    with pytest.raises(OSError, match="missing or modified"):
        manager.load_result(record.id)


def test_from_run_record_uses_persisted_identity(monkeypatch, tmp_path):
    _manager, record, output = _stored_run(tmp_path)

    class Reader:
        def __init__(self, parameters, calculator, directory, task_name=None):
            assert parameters is None
            assert task_name == "DOS"
            assert str(directory) == str(tmp_path)

        def read_from_file(self, filename, error=None, return_code=0, input_file=None):
            assert str(filename) == str(output)
            assert return_code == 0
            assert input_file is None
            return type("Result", (), {})()

    monkeypatch.setattr(reader_registry, "class_for_task_name", lambda name: Reader)
    result = TaskResult.from_run_record(record)
    assert result.run_record is record


def test_from_run_delegates_to_explicit_manager():
    expected = object()

    class Manager:
        def load_result(self, run):
            assert run == "run-id"
            return expected

    assert TaskResult.from_run("run-id", Manager()) is expected


def test_from_run_uses_default_project_store(monkeypatch, tmp_path):
    from ase2sprkkr.execution import project_manager

    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Test")
    manager = ExecutionManager(project_repository=projects, project=project)
    assert manager.store.directory == tmp_path / "projects_runs"

    output = tmp_path / "saved.out"
    output.write_text("saved result")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    manager.store.save(record)
    manager.complete(
        record, state=RunState.FINISHED, return_code=0, artifacts_available=True
    )

    monkeypatch.setattr(project_manager, "ProjectRepository", lambda: projects)
    seen = []

    def parse(cls, persisted):
        seen.append(persisted.id)
        return persisted

    monkeypatch.setattr(TaskResult, "from_run_record", classmethod(parse))
    restored = TaskResult.from_run(record.id)
    assert restored.id == record.id
    assert seen == [record.id]


def test_attach_uses_task_result_reconstruction(monkeypatch, tmp_path):
    manager, record, output = _stored_run(tmp_path)
    captured = {}
    runner = object()

    def attach(run, *, process_result, callback=None, output_listener=None):
        captured["run"] = run
        captured["process_result"] = process_result
        return runner

    monkeypatch.setattr(manager, "attach", attach)
    seen = []

    def parse(cls, persisted, **kwargs):
        seen.append((persisted, kwargs))
        return type("Result", (), {})()

    monkeypatch.setattr(TaskResult, "from_run_record", classmethod(parse))
    assert manager.attach_sprkkr(record.id) is runner

    from ase2sprkkr.execution.base import ExecutionOutcome

    outcome = ExecutionOutcome(
        input_file=None, output_file=str(output), error_file=None, return_code=0
    )
    result = captured["process_result"](outcome)
    assert result.run_record.id == record.id
    assert len(seen) == 1
    assert seen[0][0].id == record.id
    assert seen[0][1] == {
        "output_file": str(output), "error_file": None, "return_code": 0
    }
