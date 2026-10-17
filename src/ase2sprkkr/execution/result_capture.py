"""Parser-free output sink used while persistent execution is finalizing."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass


@dataclass
class CapturedResult:
    """Minimal execution result used when reconnecting without a task reader."""

    output_file: object = None
    input_file: object = None
    return_code: int | None = None

    def complete(self, _potential, return_code):
        self.return_code = return_code
        return self


class PersistentOutputCapture:
    """Consume persistent output without validating its contents.

    A live calculation may retain its original task reader for the later check.
    A reattached calculation deliberately does not need one: finalization only
    captures/retrieves files, and the task reader is reconstructed after the run
    has reached AVAILABLE.
    """

    def __init__(self, reader=None):
        self.reader = reader

    @property
    def print_output(self):
        return False if self.reader is None else self.reader.print_output

    def create_result(self, output_file, input_file=None):
        if self.reader is None:
            return CapturedResult(output_file, input_file)
        return self.reader.create_result(output_file, input_file)

    @staticmethod
    async def _drain(stream):
        while await stream.read(65536):
            pass
        return None

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(self._drain(stdout), self._drain(stderr))

    @staticmethod
    def complete_result(result, output, error, return_code):
        result.complete(None, return_code)
        return result
