""" Density of states (DOS) reader and result. """

from ..task_result import TaskResult, KkrProcessRunner, OutputFileResultValue
from .default import DefaultOutputReader
from ...common.decorators import cached_property
from ...output_files.output_files import OutputFile
import os


class DosResult(TaskResult):
  """ Density of states (DOS) results provides access to the computed
  density of states in the DOS output file, using the :py:attr:dos property. """

  @cached_property
  def dos_filename(self):
      """ Density of states file name """
      fname = self.input_parameters.CONTROL.DATASET() + '_DOS.dos'
      if self.directory:
         fname = os.path.join(self.directory, fname)
      return fname

  @cached_property
  def dos(self):
      """ The computed density of states. """
      return self.output_values['dos']()

  @cached_property
  def output_values(self):
    return {
        'dos': OutputFileResultValue('Density of states', 'dos', self.dos_filename)
        }



class DosProcessRunner(KkrProcessRunner):
  """ ARPES task output reader currently do nothing, just have a special
  result, that allow easy acces to spc output file """

  result_class = DosResult
  reader_class = DefaultOutputReader
