"""Result restoration stays in ASE2SPRKKR; the generic manager never parses."""

from types import SimpleNamespace

import pytest
from kkrunner import ExecutionOutcome, FileEntry, FileManifest, RunState, RunStore
from kkrunner.runs.manager import RunManager
from ase2sprkkr.outputs import project_results
from ase2sprkkr.outputs.task_result import TaskResult


def stored_run(tmp_path, state=RunState.FINISHED):
    manager = RunManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "saved.out"
    output.write_text("saved output")
    record = manager.new_record(
        directory=tmp_path,
        metadata={"ase2sprkkr": {"task": "DOS"}},
        processing_required=True,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.complete(record, state=state, return_code=0, artifacts_available=True)
    return manager, record, output


def test_generic_manager_has_no_application_methods(tmp_path):
    manager, _, _ = stored_run(tmp_path)
    assert not hasattr(manager, "load_result")
    assert not hasattr(manager, "attach_sprkkr")


def test_live_record_and_history_use_same_reader_factory(monkeypatch, tmp_path):
    manager, record, output = stored_run(tmp_path)
    calls = []

    class Reader:
        input_parameters = None

        def read_from_file(self, filename, error=None, **kwargs):
            calls.append(str(filename))
            return SimpleNamespace()

    def reader_factory(owner, current):
        assert current.id == record.id
        return Reader()

    monkeypatch.setattr(project_results, "_reader_for_record", reader_factory)
    first = TaskResult.from_run_record(record)
    second = project_results.load_project_result(manager, record.id)
    assert calls == [str(output), str(output)]
    assert first.run_record == second.run_record == record


def test_history_rejects_modified_output(tmp_path):
    manager, record, output = stored_run(tmp_path)
    output.write_text("overwritten")
    with pytest.raises(project_results.ProjectResultUnavailableError, match="overwritten"):
        project_results.load_project_result(manager, record.id)


def test_from_run_uses_application_adapter(monkeypatch):
    expected, manager = object(), object()
    monkeypatch.setattr(
        project_results,
        "load_project_result",
        lambda owner, run: expected if owner is manager and run == "run-id" else None,
    )
    assert TaskResult.from_run("run-id", manager) is expected


def test_attach_builds_processor_in_application_layer(monkeypatch, tmp_path):
    manager, record, output = stored_run(tmp_path, state=RunState.POSTPROCESSING)
    expected = object()
    monkeypatch.setattr(
        project_results,
        "_reader_for_record",
        lambda owner, run: SimpleNamespace(read_from_file=lambda *args, **kwargs: expected),
    )
    captured = {}

    def attach(run, **kwargs):
        captured.update(kwargs)
        return "attached"

    monkeypatch.setattr(manager, "attach", attach)
    assert project_results.attach_project_run(manager, record.id) == "attached"
    assert captured["process_result"](ExecutionOutcome(output_file=str(output), return_code=0)) is expected
