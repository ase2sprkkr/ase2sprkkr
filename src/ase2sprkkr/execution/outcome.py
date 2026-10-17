"""Parser-free execution results used below the task-result layer."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
from typing import Optional


@dataclass
class ExecutionOutcome:
    """Facts produced by execution before task-specific postprocessing.

    Execution backends own process/job lifecycle and local artifact availability.
    They deliberately do not construct a ``TaskResult``.  The higher-level
    postprocessing layer turns this outcome plus the local artifacts into the
    task-specific result.
    """

    output_file: object = None
    input_file: object = None
    return_code: Optional[int] = None


class ExecutionCapture:
    """Drain process streams without interpreting SPR-KKR output."""

    def __init__(self, *, print_output=False):
        self.print_output = bool(print_output)

    @staticmethod
    def create_result(output_file, input_file=None):
        return ExecutionOutcome(output_file=output_file, input_file=input_file)

    @staticmethod
    async def _drain(stream):
        while await stream.read(65536):
            pass
        return None

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(self._drain(stdout), self._drain(stderr))

    @staticmethod
    def complete_result(result, _output, _error, return_code):
        result.return_code = return_code
        return result
