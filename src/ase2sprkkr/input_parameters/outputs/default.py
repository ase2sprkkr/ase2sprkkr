from ..task_result import TaskResult, OutputReader

class DefaultResult(TaskResult):
      pass

class DefaultOutputReader(OutputReader):

  async def read_output(self, stdout):

      #just read the whole output
      while await stdout.readline():
        pass

  result_class = DefaultResult

