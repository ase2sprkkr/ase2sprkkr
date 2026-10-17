"""Parser-free output sink used while persistent execution is finalizing."""

from __future__ import annotations

import asyncio


class PersistentOutputCapture:
    """Adapt a task reader to execution without validating its output.

    Persistent backends still need a result placeholder and a consumer for the
    replayed stdout/stderr streams. Parsing those streams here would make parser
    validity part of execution finalization. This adapter therefore only drains
    them; the original reader is retained by :class:`PersistentRunner` and used
    later for the explicit AVAILABLE -> CHECKED phase.
    """

    def __init__(self, reader):
        self.reader = reader

    @property
    def print_output(self):
        return self.reader.print_output

    def create_result(self, output_file, input_file=None):
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
