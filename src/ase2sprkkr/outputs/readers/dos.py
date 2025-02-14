""" Density of states (DOS) reader and result. """

from ..task_result import TaskResult, KkrProcess
from .default import DefaultOutputReader
from ...common.decorators import cached_property
from ...output_files.output_files import OutputFile
import os


class DosResult(TaskResult):
  """ Density of states (DOS) results provides access to the computed
  density of states in the DOS output file, using the :py:attr:dos property. """

  @cached_property
  def dos_filename(self):
      """ New (output) potential file name """
      fname = self.input_parameters.CONTROL.DATASET() + '_DOS.dos'
      if self.directory:
         fname = os.path.join(self.directory, fname)
      return fname

  @cached_property
  def dos(self):
      """ The new (output) potential - that contains the converged charge density etc. """
      return OutputFile.from_file(self.dos_filename, try_only='dos')


class DosProcess(KkrProcess):
  """ ARPES task output reader currently do nothing, just have a special
  result, that allow easy acces to spc output file """

  result_class = DosResult
  reader_class = DefaultOutputReader
