""" This module contains classes, used by parsers of the output files """

from ..common.decorators import cached_property, add_to_signature
from ..common.process_output_reader import BaseProcessOutputReader
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


class OutputReader(BaseProcessOutputReader):

  """ Process reader, that construct (a descendant of) InputParametersResult as a result.
      Subclasses should specify result_class class property.
  """
  @add_to_signature(BaseProcessOutputReader.__init__)
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

  def result(self, *args):
      return self.result_class(self.input_parameters, self.calculator, self.directory, *args)
