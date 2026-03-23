""" The ARPES reader and result. """

from ..task_result import TaskResult, KkrProcessRunner, OutputFileResultValue
from .default import DefaultOutputReader
from ...common.decorators import cached_property
import os


class ArpesResult(TaskResult):
  """ ARPES result provides acces to the generated SPC output file
  containing the results of Angle resolved photoelectron electroscopy,
  using the :py:attr:`~spc` property. """

  @cached_property
  def spc_filename(self):
      """ Spectroscopy results file name """
      fname = self.input_parameters.CONTROL.DATASET() + '_ARPES_data.spc'
      if self.directory:
         fname = os.path.join(self.directory, fname)
      return fname

  @cached_property
  def spc(self):
      """ Spectroscopy results. """
      return self.output_values['spc']()

  @property
  def output_values(self):
    return {
        'spc': OutputFileResultValue('Spectroscopy', 'spc', self.spc_filename)
        }


class ArpesProcessRunner(KkrProcessRunner):
  """ ARPES task have no special reader, it just have a special
  result, that allow easy acces to spc output file """

  result_class = ArpesResult
  reader_class = DefaultOutputReader
