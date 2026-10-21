"""Terminal lifecycle semantics are consistent for all execution backends."""

import pytest

from ase2sprkkr.execution.base import RunState


@pytest.mark.parametrize("state", tuple(RunState))
def test_terminal_state_has_single_definition(state):
    assert state.is_terminal is (state in {
        RunState.FINISHED, RunState.FAILED, RunState.CANCELLED,
    })


def test_unknown_is_recoverable_and_never_terminal():
    assert not RunState.UNKNOWN.is_terminal
    assert not RunState.POSTPROCESSING.is_terminal
