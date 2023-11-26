from ..task_result import TaskResult, KkrProcess
from ...common.process_output_reader import ProcessOutputReader


class DefaultResult(TaskResult):
      pass


class DefaultOutputReader(ProcessOutputReader):

  async def read_output(self, stdout):

      # just read the whole output
      while await stdout.readline():
        pass


class DefaultProcess(KkrProcess):

  result_class = DefaultResult
  reader_class = DefaultOutputReader
