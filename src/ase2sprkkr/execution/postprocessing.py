"""Generic client-side result processing layered over raw execution."""

from __future__ import annotations

import asyncio
from threading import Event, Lock, Thread

from .base import FailurePhase, RunState, SubmissionNotStarted


class PostprocessingRunner:
    """Convert a raw execution outcome into the caller-visible result.

    The wrapped runner owns submission, execution and artifact availability.
    This layer owns only the ``POSTPROCESSING`` phase and invokes one supplied
    task-agnostic callable after successful raw execution.
    """

    def __init__(self, runner, process_result, callback=None):
        self._runner = runner
        self._process_result = process_result
        self.callback = callback
        self.result = None

        self._state = RunState.NEW
        self._failure_phase = None
        self._failure_message = None
        self._exception = None
        self._callback_exception = None
        self._started = False
        self._processing = False
        self._settled = False
        self._done = Event()
        self._lock = Lock()
        self._thread = None

    @property
    def spec(self):
        """Expose the immutable launch specification before submission."""
        return self._runner.spec

    @property
    def executor(self):
        """Expose the raw executor without leaking task processing into it."""
        return self._runner.executor

    @property
    def backend_id(self):
        return self._runner.backend_id

    @property
    def execution(self):
        return getattr(self._runner, 'execution', None)

    @property
    def files(self):
        return getattr(self._runner, 'files', None)

    @property
    def exception(self):
        return self._exception or self._runner.exception

    @property
    def failure_phase(self):
        return self._failure_phase or getattr(self._runner, "failure_phase", None)

    @property
    def failure_message(self):
        return self._failure_message or getattr(self._runner, "failure_message", None)

    @property
    def return_code(self):
        value = getattr(self.result, "return_code", None)
        if value is not None:
            return value
        return getattr(self._runner, "return_code", None)

    @property
    def state(self):
        if self._settled:
            return self._state
        if self._processing:
            return RunState.POSTPROCESSING
        state = self._runner.state
        if state == RunState.FINISHED and self._started:
            # Raw execution is complete, but the public lifecycle is not.
            return RunState.POSTPROCESSING
        return state

    def add_output_listener(self, callback):
        return self._runner.add_output_listener(callback)

    def add_output_error_listener(self, callback):
        return self._runner.add_output_error_listener(callback)

    def _finish(self):
        try:
            outcome = self._runner.wait()
            raw_state = self._runner.state
            if raw_state == RunState.CANCELLED:
                self.result = outcome
                self._state = RunState.CANCELLED
                return
            if raw_state == RunState.FAILED:
                self.result = outcome
                self._state = RunState.FAILED
                self._failure_phase = (
                    getattr(self._runner, "failure_phase", None)
                    or FailurePhase.EXECUTION
                )
                self._failure_message = getattr(
                    self._runner, "failure_message", None
                )
                return
            if raw_state == RunState.UNKNOWN:
                self.result = outcome
                self._state = RunState.UNKNOWN
                # An unresolved backend may still be executing. It is
                # recoverable and must not finalize as a failed project run.
                self._failure_message = getattr(
                    self._runner, "failure_message", None
                )
                return
            if raw_state != RunState.FINISHED:
                raise RuntimeError(
                    f"Raw runner completed wait() in unexpected state {raw_state.value}"
                )

            self._processing = True
            try:
                self.result = self._process_result(outcome)
            except BaseException as exc:
                self._exception = exc
                self._failure_phase = FailurePhase.POSTPROCESSING
                self._failure_message = str(exc)
                self._state = RunState.FAILED
                return
            finally:
                self._processing = False
            self._state = RunState.FINISHED
        except BaseException as exc:
            self._exception = exc
            self._failure_phase = (
                getattr(self._runner, "failure_phase", None)
                or FailurePhase.EXECUTION
            )
            self._failure_message = (
                getattr(self._runner, "failure_message", None) or str(exc)
            )
            # A failed wait() can itself mean that the backend is unreachable.
            # Never replace the original error by a second status-query error.
            try:
                cancelled = self._runner.state == RunState.CANCELLED
            except Exception:
                cancelled = False
            self._state = RunState.CANCELLED if cancelled else RunState.FAILED
        finally:
            # Make terminal status visible during the persistence callback.
            self._settled = True
            try:
                if self.callback is not None:
                    self.callback(self.result)
            except BaseException as exc:
                self._callback_exception = exc
            self._done.set()

    def _start_completion(self):
        self._started = True
        self._thread = Thread(target=self._finish, daemon=True)
        self._thread.start()

    def prepare_submission(self):
        """Delegate side-effect-free-of-process launch preparation to raw runner."""
        prepare = getattr(self._runner, "prepare_submission", None)
        if prepare is None:
            return False
        prepare()
        return True

    def submit(self):
        with self._lock:
            if self._started or self._runner.state != RunState.NEW:
                raise RuntimeError(f"Run is already {self.state.value}")
            try:
                self._runner.submit()
            except SubmissionNotStarted:
                raise
            self._start_completion()

    def start(self):
        """Start a new run or monitoring of an attached raw runner."""
        with self._lock:
            if self._started:
                return
            if self._runner.state == RunState.NEW:
                try:
                    self._runner.submit()
                except SubmissionNotStarted:
                    raise
            else:
                start = getattr(self._runner, "start", None)
                if start is not None:
                    start()
            self._start_completion()

    def wait(self):
        if not self._started:
            raise RuntimeError("Run has not been started")
        self._done.wait()
        if self._exception is not None:
            raise self._exception
        if self._callback_exception is not None:
            raise self._callback_exception
        return self.result

    def run(self):
        if not self._started:
            self.start()
        return self.wait()

    async def run_async(self):
        if not self._started:
            self.start()
        return await asyncio.to_thread(self.wait)

    def cancel(self):
        return self._runner.cancel()
