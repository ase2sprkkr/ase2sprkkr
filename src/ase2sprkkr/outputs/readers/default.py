from ..task_result import TaskResult, KkrProcess
from ..sprkkr_output_reader import SprKkrOutputReader


class DefaultResult(TaskResult):
      pass


class DefaultOutputReader(SprKkrOutputReader):

  async def read_output(self, stdout, result):
      await self.read_commons(stdout, result)
      # just consume the whole rest of output
      while await stdout.readline():
        pass


class DefaultProcess(KkrProcess):

  result_class = DefaultResult
  reader_class = DefaultOutputReader
