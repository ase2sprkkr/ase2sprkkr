import asyncio
import os
import sys
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution.base import (
    ExecutionSpec,
    RunState,
    SubmissionNotStarted,
)
from ase2sprkkr.execution.detached import _terminate_process_group
from ase2sprkkr.execution.endpoint import endpoint_for_executor
from ase2sprkkr.execution.local import LocalExecutor, LocalRunner
from ase2sprkkr.execution.manager import (
    ExecutionManager,
    FileRetrievalState,
    PersistentRunner,
    ResultFilesRecord,
    RunStore,
    SubmissionKey,
)
from ase2sprkkr.execution.transport import AtomicLocalFile, LocalTransport


class _Reader:
    print_output = False

    def __init__(self, parser_error=None):
        self.parser_error = parser_error

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None)

    async def parse_streams(self, stdout, stderr, result):
        output, error = await asyncio.gather(stdout.read(), stderr.read())
        if self.parser_error is not None:
            raise self.parser_error
        return output, error

    def complete_result(self, result, output, error, return_code):
        result.return_code = return_code
        return result


def _local_spec(tmp_path, code="print('ok')"):
    return ExecutionSpec(
        argv=(sys.executable, "-c", code),
        input_file=None,
        output_file=(tmp_path / "out.txt").open("w+b"),
        error_file=(tmp_path / "err.txt").open("w+b"),
        directory=str(tmp_path),
    )


def test_output_observer_failure_is_reported_but_does_not_fail_process(tmp_path):
    runner = LocalRunner(_Reader(), _local_spec(tmp_path))
    observer_errors = []

    def broken_observer(*_):
        raise ValueError("observer failed")

    runner.add_output_listener(broken_observer)
    runner.add_output_error_listener(
        lambda callback, exc: observer_errors.append((callback, exc))
    )

    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert len(observer_errors) == 1
    assert observer_errors[0][0] is broken_observer
    assert str(observer_errors[0][1]) == "observer failed"


def test_parser_failure_does_not_change_local_execution_state(tmp_path):
    runner = LocalRunner(
        _Reader(ValueError("parser failed")),
        _local_spec(tmp_path),
    )

    runner.submit()
    with pytest.raises(ValueError, match="parser failed"):
        runner.wait()

    assert runner.state == RunState.FINISHED


def test_callback_failure_does_not_change_local_execution_state(tmp_path):
    runner = LocalRunner(
        _Reader(),
        _local_spec(tmp_path),
        callback=lambda result: (_ for _ in ()).throw(ValueError("callback failed")),
    )

    runner.submit()
    with pytest.raises(ValueError, match="callback failed"):
        runner.wait()

    assert runner.state == RunState.FINISHED


class _NotStartedRunner:
    result = SimpleNamespace(return_code=None)
    reference = None
    state = RunState.NEW

    def submit(self):
        raise SubmissionNotStarted("input staging failed")

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _UncertainRunner(_NotStartedRunner):
    def submit(self):
        raise OSError("connection lost during submit")


def _persistent_wrapper(tmp_path, runner):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = LocalExecutor()
    endpoint = endpoint_for_executor(executor)
    submission = SubmissionKey(endpoint, str(tmp_path), "SCF")
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        endpoint=endpoint,
        submission=submission,
    )
    return manager, PersistentRunner(manager, executor, runner, record)


def test_definite_pre_submit_failure_becomes_execution_failed(tmp_path):
    manager, runner = _persistent_wrapper(tmp_path, _NotStartedRunner())

    with pytest.raises(SubmissionNotStarted, match="input staging failed"):
        runner.submit()

    record = manager.store.get(runner.id)
    assert record.state == RunState.EXECUTION_FAILED
    assert not hasattr(record, "backend_state")
    assert record.reference is None
    assert record.message == "input staging failed"


def test_uncertain_submit_without_reference_stays_submitting(tmp_path):
    manager, runner = _persistent_wrapper(tmp_path, _UncertainRunner())

    with pytest.raises(OSError, match="connection lost during submit"):
        runner.submit()

    record = manager.store.get(runner.id)
    assert record.state == RunState.SUBMITTING
    assert record.reference is None


def test_file_retrieval_failure_does_not_change_finished_state(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        result_files=ResultFilesRecord(output_file=str(tmp_path / "out")),
    )
    record = manager.store.save(record)
    record = manager.transition(record, state=RunState.FINISHED)

    record = manager.update_file_retrieval(
        record, FileRetrievalState.FAILED, "ssh disconnected"
    )

    assert record.state == RunState.FINISHED
    assert record.result_files.retrieval_state == FileRetrievalState.FAILED
    assert record.result_files.retrieval_message == "ssh disconnected"


def test_atomic_local_file_preserves_old_destination_until_commit(tmp_path):
    destination = tmp_path / "output"
    destination.write_bytes(b"old")
    replay = AtomicLocalFile(destination)
    replay.write(b"new")
    replay.flush()

    assert destination.read_bytes() == b"old"
    replay.discard()
    assert destination.read_bytes() == b"old"

    replay = AtomicLocalFile(destination)
    replay.write(b"new")
    replay.commit()
    assert destination.read_bytes() == b"new"


def test_local_transport_failed_copy_preserves_destination(tmp_path):
    destination = tmp_path / "output"
    destination.write_bytes(b"old")

    with pytest.raises(FileNotFoundError):
        asyncio.run(LocalTransport().put_file(tmp_path / "missing", destination))

    assert destination.read_bytes() == b"old"


def test_process_group_cancel_checks_group_not_only_leader(tmp_path, monkeypatch):
    calls = []

    def killpg(pid, sig):
        calls.append((pid, sig))
        if sig == 0:
            raise ProcessLookupError

    monkeypatch.setattr(os, "killpg", killpg)
    monkeypatch.setattr(
        os,
        "kill",
        lambda *args: (_ for _ in ()).throw(AssertionError("leader-only check used")),
    )
    monkeypatch.setattr(
        "ase2sprkkr.execution.detached._read_state",
        lambda directory, run_id: {"state": "cancelled"},
    )

    state = _terminate_process_group(tmp_path, "run", 1234, timeout=1.0)

    assert state["state"] == "cancelled"
    assert calls[0][0] == 1234
    assert any(sig == 0 for _, sig in calls)
