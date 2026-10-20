"""ProjectRunner state observation and exactly-once finalization."""

from types import SimpleNamespace

from ..base import RunState
from ..manager import ProjectRunner


class _Manager:
    def __init__(self):
        self.completions = 0
        self.transitions = 0

    def complete(self, record, **kwargs):
        self.completions += 1
        return SimpleNamespace(id=record.id, state=kwargs["state"])

    def transition(self, record, **kwargs):
        self.transitions += 1
        return record


class _Runner:
    state = RunState.RUNNING
    backend_id = "backend-id"
    result = SimpleNamespace(return_code=0)
    exception = None
    failure_phase = None
    failure_message = None
    return_code = 0
    files = None
    execution = None

    def wait(self):
        self.state = RunState.FINISHED
        return self.result


class _Executor:
    def __init__(self):
        self.acknowledgements = 0

    def acknowledge(self, record):
        self.acknowledgements += 1


def test_reading_state_does_not_write_to_store():
    manager, raw, executor = _Manager(), _Runner(), _Executor()
    record = SimpleNamespace(id="run", state=RunState.RUNNING,
                             updated_at="2026-10-20T10:00:00+00:00",
                             execution=SimpleNamespace(submitted=lambda x: x))
    run = ProjectRunner(manager, executor, raw, record)
    assert run.state == RunState.RUNNING
    assert manager.transitions == 0


def test_finalization_is_once_even_if_wait_and_callback_both_request_it():
    manager, raw, executor = _Manager(), _Runner(), _Executor()
    record = SimpleNamespace(id="run", state=RunState.RUNNING,
                             updated_at="2026-10-20T10:00:00+00:00",
                             execution=SimpleNamespace(submitted=lambda x: x))
    run = ProjectRunner(manager, executor, raw, record)
    raw.state = RunState.FINISHED
    run._sync_completion(raw.result)
    run._sync_completion(raw.result)
    run.wait()
    assert manager.completions == 1
    assert executor.acknowledgements == 1


def test_unknown_observation_is_written_once_but_not_finalized():
    manager, raw, executor = _Manager(), _Runner(), _Executor()
    record = SimpleNamespace(id="run", state=RunState.RUNNING,
                             updated_at="2026-10-20T10:00:00+00:00",
                             execution=SimpleNamespace(submitted=lambda x: x))
    run = ProjectRunner(manager, executor, raw, record)
    raw.state = RunState.UNKNOWN
    raw.failure_message = "scheduler disconnected"

    run._sync_completion(raw.result)
    run._sync_completion(raw.result)
    run.wait = lambda: run._sync_completion(raw.result)
    run.wait()

    assert manager.transitions == 1
    assert manager.completions == 0
    assert executor.acknowledgements == 0


def test_unknown_observation_does_not_fire_completion_callback():
    from ..manager import ExecutionManager

    manager, raw, executor = _Manager(), _Runner(), _Executor()
    record = SimpleNamespace(id="run", state=RunState.RUNNING,
                             updated_at="2026-10-20T10:00:00+00:00",
                             execution=SimpleNamespace(submitted=lambda x: x))
    run = ProjectRunner(manager, executor, raw, record)
    delivered = []
    raw.state = RunState.UNKNOWN
    ExecutionManager._complete(run, raw.result, delivered.append)
    assert delivered == [raw.result]

    raw.state = RunState.FINISHED
    new_run = ProjectRunner(manager, executor, raw, record)
    ExecutionManager._complete(new_run, raw.result, delivered.append)
    assert delivered == [raw.result, raw.result]

