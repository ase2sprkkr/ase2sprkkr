"""Common parent for all specialized readers and outputs and
default reader for the tasks without specialized reader and output."""

from ..task_result import TaskResult, KkrOutputReader
from ..sprkkr_output_reader import SprKkrOutputParser


class DefaultResult(TaskResult):
    pass


class DefaultOutputParser(SprKkrOutputParser):
    async def read_output(self, stdout, result):
        await self.read_commons(stdout, result)
        # just consume the whole rest of output

        result.output_lines = []
        async for line in stdout:
            result.output_lines.append(line.decode("utf8").rstrip())
        return result


class DefaultOutputReader(KkrOutputReader):
    result_class = DefaultResult
    parser_class = DefaultOutputParser
