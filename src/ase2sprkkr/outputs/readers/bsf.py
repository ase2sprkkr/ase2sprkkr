""" The Bloch spectral functions (BSF) reader and result."""

from ..task_result import TaskResult, KkrProcessRunner, OutputFileResultValue
from .default import DefaultOutputReader
from ...common.decorators import cached_property
from ...output_files.output_files import OutputFile


class BsfResult(TaskResult):
  """ BSF result provides access to the computed Bloch spectral functions in the
  BSF output file using the :py:attr:`~bsf` property."""

  @cached_property
  def bsf_filename(self):
      """ New (output) potential file name """
      return self.path_to('Bloch-SF')

  @cached_property
  def bsf(self):
      """ The new (output) potential - that contains the converged charge density etc. """
      self.output_values['bsf'].output_file()

  @cached_property
  def output_values(self):
    return {
        'bsf': OutputFileResultValue('Bloch spectral function', 'bsf', self.bsf_filename)
        }

class BsfProcessRunner(KkrProcessRunner):
  """ ARPES task output reader currently do nothing, just have a special
  result, that allow easy acces to BSF output file """

  result_class = BsfResult
  reader_class = DefaultOutputReader
