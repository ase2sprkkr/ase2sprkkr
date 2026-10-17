"""Raw-execution facade over the durable run-history implementation."""

from __future__ import annotations

from .manager import ExecutionManager as _ExecutionManager
from .postprocessing import PostprocessingRunner


class ExecutionManager(_ExecutionManager):
    """Persist a raw execution runner plus client-side postprocessing.

    The durable store remains task-agnostic.  ``postprocess`` is a caller-owned
    callable that turns an execution outcome into the public result after local
    artifacts are available; it is never passed to an execution backend.
    """

    def create_run(
        self,
        executor,
        spec,
        *,
        postprocess,
        print_output=False,
        task,
        input_definition=None,
        dataset=None,
        callback=None,
        run_id=None,
    ):
        endpoint = __import__(
            "ase2sprkkr.execution.endpoint",
            fromlist=["endpoint_for_executor"],
        ).endpoint_for_executor(executor)
        from .manager import SubmissionKey

        submission = SubmissionKey(
            endpoint=endpoint,
            directory=executor.submission_directory(spec),
            task=task,
            dataset=dataset,
        )
        record = self.new_record(
            task=task,
            input_definition=input_definition,
            dataset=dataset,
            directory=spec.directory,
            endpoint=endpoint,
            submission=submission,
            result_files=self._result_files_for_spec(spec),
            run_id=run_id,
        )
        self.register_executor(executor, endpoint)

        def factory(completed):
            raw = executor.create_persistent_run(
                spec,
                run_id=record.id,
                print_output=print_output,
            )
            return PostprocessingRunner(raw, postprocess, callback=completed)

        return self._wrap_runner(executor, record, factory, callback)

    def attach(
        self,
        run,
        *,
        postprocess,
        print_output=False,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        """Reconnect raw execution first; reconstruct the result client-side."""
        from .endpoint import endpoint_for_executor

        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        if executor is None:
            executor = self.executor_for(record)
        else:
            endpoint = record.endpoint or endpoint_for_executor(executor)
            self.register_executor(executor, endpoint, record.reference)

        def factory(completed):
            raw = executor.attach(
                record.reference,
                output_listener=output_listener,
                start=False,
                print_output=print_output,
            )
            wrapped = PostprocessingRunner(raw, postprocess, callback=completed)
            # Start through the wrapper so POSTPROCESSING remains observable and
            # the result callback is tied to the complete client-side lifecycle.
            wrapped.submit()
            return wrapped

        return self._wrap_runner(executor, record, factory, callback)
