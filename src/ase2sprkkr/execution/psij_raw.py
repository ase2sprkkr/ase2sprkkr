"""Parser-free public PSI/J executor built on the scheduler transport runner."""

from __future__ import annotations

from .outcome import ExecutionCapture
from .psij import PsijExecutor as _PsijExecutor
from .psij import PsijRunner


class PsijExecutor(_PsijExecutor):
    """PSI/J executor whose public boundary contains no task output reader."""

    def create_run(self, spec, callback=None, *, print_output=False):
        return PsijRunner(
            self,
            ExecutionCapture(print_output=print_output),
            spec,
            callback=callback,
        )

    def create_persistent_run(
        self,
        spec,
        *,
        run_id,
        callback=None,
        print_output=False,
    ):
        return PsijRunner(
            self,
            ExecutionCapture(print_output=print_output),
            spec,
            callback=callback,
            run_id=run_id,
        )

    def attach(
        self,
        reference,
        callback=None,
        output_listener=None,
        start=True,
        *,
        print_output=False,
    ):
        return super().attach(
            reference,
            ExecutionCapture(print_output=print_output),
            callback=callback,
            output_listener=output_listener,
            start=start,
        )
