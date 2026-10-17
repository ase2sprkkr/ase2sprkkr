"""Raw-execution facade over the durable run-history implementation."""

from __future__ import annotations

from .endpoint import endpoint_for_executor
from .manager import ExecutionManager as _ExecutionManager
from .manager import SubmissionKey
from .postprocessing import PostprocessingRunner, reader_postprocessor


class ExecutionManager(_ExecutionManager):
    """Persist a raw execution runner plus client-side postprocessing.

    New callers pass ``spec`` plus ``postprocess=...``.  The historical
    ``(reader, spec)`` form remains temporarily accepted at this facade only;
    the reader is immediately converted to a client-side postprocessor and is
    never forwarded to the executor.
    """

    @staticmethod
    def _normalize_postprocess(reader_or_spec, spec, postprocess):
        if spec is None:
            return reader_or_spec, postprocess
        if postprocess is not None:
            raise TypeError("postprocess cannot be combined with legacy reader argument")
        return spec, reader_postprocessor(reader_or_spec, spec)

    def create_run(
        self,
        executor,
        reader_or_spec,
        spec=None,
        *,
        postprocess=None,
        print_output=False,
        task,
        input_definition=None,
        dataset=None,
        callback=None,
        run_id=None,
    ):
        spec, postprocess = self._normalize_postprocess(
            reader_or_spec, spec, postprocess
        )
        if postprocess is None:
            raise TypeError("postprocess is required for persistent execution")

        endpoint = endpoint_for_executor(executor)
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
        reader=None,
        *,
        postprocess=None,
        print_output=False,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        """Reconnect raw execution first; reconstruct the result client-side."""
        record = self.store.get(run) if isinstance(run, str) else run
        if record.reference is None:
            raise RuntimeError(f"Run {record.id} has no backend reference")
        if executor is None:
            executor = self.executor_for(record)
        else:
            endpoint = record.endpoint or endpoint_for_executor(executor)
            self.register_executor(executor, endpoint, record.reference)

        if postprocess is None:
            if reader is None:
                raise TypeError("postprocess or legacy reader is required")
            # Reattached raw runners reconstruct their ExecutionSpec from the
            # persisted reference, so build an equivalent lightweight spec here.
            from .base import ExecutionSpec

            spec = ExecutionSpec(
                argv=(),
                input_file=record.reference.input_file,
                output_file=record.reference.output_file,
                error_file=record.reference.error_file,
                directory=record.directory,
            )
            postprocess = reader_postprocessor(reader, spec)

        def factory(completed):
            raw = executor.attach(
                record.reference,
                output_listener=output_listener,
                start=False,
                print_output=print_output,
            )
            wrapped = PostprocessingRunner(raw, postprocess, callback=completed)
            wrapped.submit()
            return wrapped

        return self._wrap_runner(executor, record, factory, callback)
