"""Client-side task-result postprocessing layered over raw execution."""

from __future__ import annotations

import asyncio
from threading import Event, Lock, Thread

from .base import FailurePhase, RunState, SubmissionNotStarted


class PostprocessingRunner:
    """Turn a raw execution outcome into the caller-visible result.

    The wrapped runner owns execution and artifact availability.  This wrapper
    owns only task-independent orchestration around a supplied ``postprocess``
    callable.  The callable runs after successful execution, when all required
    artifacts are local, and is therefore the only layer that may construct a
    task-specific result.
    """

    def __init__(self, runner, postprocess, callback=None):
        self._runner = runner
        self._postprocess = postprocess
        self.callback = callback
        self.result = None
        self._state = RunState.NEW
        self._failure_phase = None
        self._failure_message = None
        self._exception = None
        self._callback_exception = None
        self._submitted = False
        self._processing = False
        self._done = Event()
        self._lock = Lock()
        self._thread = None

    @property
    def reference(self):
        return self._runner.reference

    @property
    def exception(self):
        return self._exception

    @property
    def failure_phase(self):
        phase = self._failure_phase
        if phase is not None:
            return phase
        return getattr(self._runner, "failure_phase", None)

    @property
    def failure_message(self):
        message = self._failure_message
        if message is not None:
            return message
        return getattr(self._runner, "failure_message", None)

    @property
    def return_code(self):
        result = self.result
        if result is not None and getattr(result, "return_code", None) is not None:
            return result.return_code
        return getattr(self._runner, "return_code", None)

    @property
    def state(self):
        if self._done.is_set():
            return self._state
        if self._processing:
            return RunState.POSTPROCESSING
        raw = self._runner.state
        if raw == RunState.FINISHED and self._submitted:
            # The execution boundary has completed, but this wrapper has not.
            return RunState.POSTPROCESSING
        return raw

    @property
    def file_retrieval_state(self):
        return getattr(self._runner, "file_retrieval_state", None)

    @property
    def file_retrieval_message(self):
        return getattr(self._runner, "file_retrieval_message", None)

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
                return

            self._processing = True
            try:
                self.result = self._postprocess(outcome)
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
            raw_phase = getattr(self._runner, "failure_phase", None)
            self._failure_phase = raw_phase or FailurePhase.EXECUTION
            self._failure_message = (
                getattr(self._runner, "failure_message", None) or str(exc)
            )
            self._state = (
                RunState.CANCELLED
                if self._runner.state == RunState.CANCELLED
                else RunState.FAILED
            )
        finally:
            try:
                if self.callback is not None:
                    self.callback(self.result)
            except BaseException as exc:
                self._callback_exception = exc
            self._done.set()

    def submit(self):
        with self._lock:
            if self._submitted or self._state != RunState.NEW:
                raise RuntimeError(f"Run is already {self.state.value}")
            try:
                self._runner.submit()
            except SubmissionNotStarted:
                raise
            self._submitted = True
            self._thread = Thread(target=self._finish, daemon=True)
            self._thread.start()

    def wait(self):
        if not self._submitted:
            raise RuntimeError("Run has not been submitted")
        self._done.wait()
        if self._exception is not None:
            raise self._exception
        if self._callback_exception is not None:
            raise self._callback_exception
        return self.result

    def run(self):
        if not self._submitted:
            self.submit()
        return self.wait()

    async def run_async(self):
        if not self._submitted:
            self.submit()
        return await asyncio.to_thread(self.wait)

    def cancel(self):
        return self._runner.cancel()
