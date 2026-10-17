"""Parser-free remote executor and explicit artifact-download lifecycle."""

from __future__ import annotations

from .base import FailurePhase, RunState
from .outcome import ExecutionCapture
from .postprocessing import PostprocessingRunner, reader_postprocessor
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
    """Remote endpoint whose backend boundary is parser-free."""

    @staticmethod
    def _capture(print_output):
        return ExecutionCapture(print_output=print_output)

    @staticmethod
    def _wrap(runner):
        return _RemoteLifecycleRunner(runner)

    def _raw_run(self, spec, callback=None, *, print_output=False, run_id=None):
        if run_id is None:
            runner = super().create_run(
                self._capture(print_output),
                spec,
                callback=callback,
            )
        else:
            runner = super().create_persistent_run(
                self._capture(print_output),
                spec,
                run_id=run_id,
                callback=callback,
            )
        return self._wrap(runner)

    def create_run(
        self,
        reader_or_spec,
        spec=None,
        callback=None,
        *,
        print_output=False,
    ):
        if spec is None:
            return self._raw_run(
                reader_or_spec,
                callback=callback,
                print_output=print_output,
            )
        reader = reader_or_spec
        raw = self._raw_run(
            spec,
            print_output=bool(getattr(reader, "print_output", print_output)),
        )
        return PostprocessingRunner(
            raw,
            reader_postprocessor(reader, spec),
            callback=callback,
        )

    def create_persistent_run(
        self,
        reader_or_spec,
        spec=None,
        *,
        run_id,
        callback=None,
        print_output=False,
    ):
        if spec is None:
            return self._raw_run(
                reader_or_spec,
                callback=callback,
                print_output=print_output,
                run_id=run_id,
            )
        reader = reader_or_spec
        raw = self._raw_run(
            spec,
            print_output=bool(getattr(reader, "print_output", print_output)),
            run_id=run_id,
        )
        return PostprocessingRunner(
            raw,
            reader_postprocessor(reader, spec),
            callback=callback,
        )

    def attach(
        self,
        reference,
        reader=None,
        callback=None,
        output_listener=None,
        start=True,
        *,
        print_output=False,
    ):
        runner = super().attach(
            reference,
            self._capture(
                bool(getattr(reader, "print_output", print_output))
            ),
            callback=None if reader is not None else callback,
            output_listener=output_listener,
            start=start if reader is None else False,
        )
        raw = self._wrap(runner)
        if reader is None:
            return raw

        from .base import ExecutionSpec

        spec = ExecutionSpec(
            argv=(),
            input_file=reference.input_file,
            output_file=reference.output_file,
            error_file=reference.error_file,
            directory=reference.directory or ".",
        )
        wrapped = PostprocessingRunner(
            raw,
            reader_postprocessor(reader, spec),
            callback=callback,
        )
        if start:
            wrapped.submit()
        return wrapped
