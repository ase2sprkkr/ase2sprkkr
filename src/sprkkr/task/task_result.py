from ..common.misc import cached_property, add_to_signature
from ..potential.potentials import Potential
from ..common.process_output_reader import BaseProcessOutputReader
import os

class TaskResult:
  """ A base class for a result of a runned task (kkrscf executable) """
  def __init__(self, task, return_code):
      self.task = task
      self.directory = task.directory
      self.return_code = return_code

  @cached_property
  def potential(self):
      fname = self.task.CONTROL.POTFIL() + '_new'
      if self.directory:
         fname = os.path.join(self.directory, fname)
      return Potential.from_file(fname)

  @cached_property
  def atoms(self):
      return self.potential.atoms


class TaskResultReader(BaseProcessOutputReader):

  """ Process reader, that construct (a descendant of) TaskResult as a result. 
      Subclasses should specify result_class class property.
  """

  @add_to_signature(BaseProcessOutputReader.__init__)
  def __init__(self, *args, task, **kwargs):
      self.task = task
      super().__init__(*args, **kwargs)

  def result(self, *args):
      return self.result_class(self.task, *args)
