from types import SimpleNamespace
from threading import Event

from ase2sprkkr.execution import RunState
from ase2sprkkr.execution.base import OutputObservers
from ase2sprkkr.execution.psij import PsijRunner


def test_psij_finalization_failure_overrides_completed_backend():
    runner = object.__new__(PsijRunner)
    runner._job = SimpleNamespace(
        status=SimpleNamespace(state=SimpleNamespace(name="COMPLETED"))
    )
    runner._exception = None
    runner._failure_phase = None
    runner._failure_message = None
    runner._return_code = None
    runner._output_observers = OutputObservers()
    runner._done = Event()
    runner.callback = None
    runner.result = object()

    async def fail():
        raise ValueError("result parsing failed")

    runner._monitor = fail
    runner._thread_main()

    assert isinstance(runner._exception, ValueError)
    assert runner.state == RunState.FAILED
    assert runner.failure_phase.value == "execution"
