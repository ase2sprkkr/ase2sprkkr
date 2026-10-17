"""Public local executor facade over the parser-free subprocess backend."""

from __future__ import annotations

from .local import LocalExecutor as _LocalExecutor
from .postprocessing import PostprocessingRunner, reader_postprocessor


class LocalExecutor(_LocalExecutor):
    """Local executor with a temporary compatibility seam for task readers."""

    def create_run(
        self,
        reader_or_spec,
        spec=None,
        callback=None,
        *,
        print_output=False,
    ):
        if spec is None:
            return super().create_run(
                reader_or_spec,
                callback=callback,
                print_output=print_output,
            )

        reader = reader_or_spec
        raw = super().create_run(
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
            return super().create_persistent_run(
                reader_or_spec,
                run_id=run_id,
                callback=callback,
                print_output=print_output,
            )

        reader = reader_or_spec
        raw = super().create_persistent_run(
            spec,
            run_id=run_id,
            print_output=bool(getattr(reader, "print_output", print_output)),
        )
        return PostprocessingRunner(
            raw,
            reader_postprocessor(reader, spec),
            callback=callback,
        )
