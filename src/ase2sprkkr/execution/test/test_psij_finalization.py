from types import SimpleNamespace
from threading import Event

from ase2sprkkr.execution import RunState
from ase2sprkkr.execution.psij import PsijRunner


def test_psij_finalization_failure_overrides_completed_backend():
    runner = object.__new__(PsijRunner)
    runner._job = SimpleNamespace(
        status=SimpleNamespace(state=SimpleNamespace(name="COMPLETED"))
    )
    runner._exception = None
    runner._finalization_failed = False
    runner._done = Event()
    runner.callback = None
    runner.result = object()

    async def fail():
        raise ValueError("result parsing failed")

    runner._monitor = fail
    runner._thread_main()

    assert isinstance(runner._exception, ValueError)
    assert runner.state == RunState.FAILED
