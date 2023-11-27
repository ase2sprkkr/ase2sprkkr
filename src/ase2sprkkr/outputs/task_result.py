""" This module contains classes, used by parsers of the output files """
from ..common.decorators import cached_property
import os


class TaskResult:
  """ A base class for a result of a runned task (kkrscf executable) """
  def __init__(self, input_parameters, calculator, directory, result, error, return_code):
      self.input_parameters = input_parameters
      self._calculator = calculator
      self.directory = directory or os.getcwd()
      self.result = result
      self.error = error
      self.return_code = return_code

  @cached_property
  def atoms(self):
      return self.potential.atoms


class KkrProcess:
  """ Class, that run a process and read its output using underlined
  process reader (see :class:`ase2sprkkr.common.process_output_reader.ProcessOutputReader`)
  and return the appropriate TaskResult.

  Descendants should define reader_class and result_class property.
  """

  def __init__(self, input_parameters, calculator, directory):
      self.input_parameters = input_parameters
      """ Input parameters, that command to read the output (thus probably the ones, that
      run the process that produced the output. It is used e.g. for determining the potential file,
      which belongs to the output.
      """
      self.calculator = calculator
      """ Calculator, that can be used for further processing of the results. """
      self.directory = directory
      """ Directory, to wich are the relative paths in the output related. """
      self.reader = self.reader_class()

  def run(self, cmd, outfile, print_output=False, directory=None, **kwargs):
      return self._wraps(*self.reader.run(cmd, outfile, print_output, directory, **kwargs))

  def _wraps(self, *args):
      return self.result_class(self.input_parameters, self.calculator, self.directory, *args)

  def read_from_file(self, output, error=None, return_code=0, print_output=False):
      return self._wraps(*self.reader.read_from_file(output, error, return_code, print_output))
