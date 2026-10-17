"""Parser-free remote executor and explicit artifact-download lifecycle."""

from __future__ import annotations

from .base import FailurePhase, RunState
from .outcome import ExecutionCapture
from .remote_endpoint import RemoteEndpointExecutor as _RemoteEndpointExecutor


class _RemoteLifecycleRunner:
    """Expose remote artifact transfer as an explicit DOWNLOADING phase."""

    def __init__(self, runner):
        self._runner = runner

    def __getattr__(self, name):
        return getattr(self._runner, name)

    @property
    def result(self):
        return self._runner.result

    @property
    def reference(self):
        return self._runner.reference

    @property
    def exception(self):
        return self._runner.exception

    @property
    def state(self):
        state = self._runner.state
        retrieval = getattr(self._runner, "file_retrieval_state", None)
        if retrieval == "failed":
            return RunState.FAILED
        if retrieval == "pending" and state in (RunState.FINISHED, RunState.FAILED):
            return RunState.DOWNLOADING
        return state

    @property
    def failure_phase(self):
        retrieval = getattr(self._runner, "file_retrieval_state", None)
        if retrieval == "failed":
            return FailurePhase.DOWNLOAD
        if self.state == RunState.FAILED:
            return FailurePhase.EXECUTION
        return None

    @property
    def failure_message(self):
        if getattr(self._runner, "file_retrieval_state", None) == "failed":
            return getattr(self._runner, "file_retrieval_message", None)
        return None


class RemoteEndpointExecutor(_RemoteEndpointExecutor):
    """Remote endpoint whose public execution API is parser-free."""

    @staticmethod
    def _capture(print_output):
        return ExecutionCapture(print_output=print_output)

    @staticmethod
    def _wrap(runner):
        return _RemoteLifecycleRunner(runner)

    def create_run(self, spec, callback=None, *, print_output=False):
        runner = super().create_run(
            self._capture(print_output),
            spec,
            callback=callback,
        )
        return self._wrap(runner)

    def create_persistent_run(
        self,
        spec,
        *,
        run_id,
        callback=None,
        print_output=False,
    ):
        runner = super().create_persistent_run(
            self._capture(print_output),
            spec,
            run_id=run_id,
            callback=callback,
        )
        return self._wrap(runner)

    def attach(
        self,
        reference,
        callback=None,
        output_listener=None,
        start=True,
        *,
        print_output=False,
    ):
        runner = super().attach(
            reference,
            self._capture(print_output),
            callback=callback,
            output_listener=output_listener,
            start=start,
        )
        return self._wrap(runner)
