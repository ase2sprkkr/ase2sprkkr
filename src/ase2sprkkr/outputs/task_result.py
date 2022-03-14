""" This module contains classes, used by parsers of the output files """

from ..common.misc import cached_property, add_to_signature
from ..common.process_output_reader import BaseProcessOutputReader

class TaskResult:
  """ A base class for a result of a runned task (kkrscf executable) """
  def __init__(self, input_parameters, calculator, result, error, return_code):
      self.input_parameters = input_parameters
      self._calculator = calculator
      self.directory = input_parameters.directory
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
  def __init__(self, input_parameters, calculator):
      self.input_parameters = input_parameters
      self.calculator = calculator

  def result(self, *args):
      return self.result_class(self.input_parameters, self.calculator, *args)
