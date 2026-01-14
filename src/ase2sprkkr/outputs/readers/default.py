""" Common parent for all specialized readers and outputs and
default reader for the tasks without specialized reader and output. """

from ..task_result import TaskResult, KkrProcessRunner
from ..sprkkr_output_reader import SprKkrOutputReader


class DefaultResult(TaskResult):
      pass


class DefaultOutputReader(SprKkrOutputReader):

  async def read_output(self, stdout, result):
      await self.read_commons(stdout, result)
      # just consume the whole rest of output

      result.output_lines = []
      async for line in stdout:
            result.output_lines.append(line.decode('utf8').rstrip())
      return result


class DefaultProcessRunner(KkrProcessRunner):

  result_class = DefaultResult
  reader_class = DefaultOutputReader
