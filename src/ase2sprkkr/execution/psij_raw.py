"""Parser-free public PSI/J executor built on the scheduler transport runner."""

from __future__ import annotations

from .outcome import ExecutionCapture
from .postprocessing import PostprocessingRunner, reader_postprocessor
from .psij import PsijExecutor as _PsijExecutor
from .psij import PsijRunner


class PsijExecutor(_PsijExecutor):
    """PSI/J executor whose backend boundary contains no task output reader."""

    def _raw_run(self, spec, callback=None, *, print_output=False, run_id=None):
        return PsijRunner(
            self,
            ExecutionCapture(print_output=print_output),
            spec,
            callback=callback,
            run_id=run_id,
        )

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
        raw = super().attach(
            reference,
            ExecutionCapture(
                print_output=bool(getattr(reader, "print_output", print_output))
            ),
            callback=None if reader is not None else callback,
            output_listener=output_listener,
            start=start if reader is None else False,
        )
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
