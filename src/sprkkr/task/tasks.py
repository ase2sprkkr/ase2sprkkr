import os
import io
import tempfile
import pkgutil
import importlib
from . import definitions
from ..common.conf_containers import RootConfContainer
from ..common.misc import lazy_value, OrderedDict

class Task(RootConfContainer):
  def __init__(self, definition, inputfile=None, outputfile=False):
      super().__init__(definition)
      self._inputfile = inputfile
      self._outputfile = inputfile if outputfile is False else outputfile

  @property
  def task_name(self):
      return self._definition.name

  def run_task_process(self, calculator, task_file, output_file, print_output=False):
      d = self._definition
      command = d.command
      if d.mpi and calculator.mpi:
           command = calculator.mpi + [ command + 'MPI' ]
      else:
           command = [ command ]

      process = d.result_reader(command, output_file, stdin = task_file, print_output=print_output)
      process.run()
      return process

  def executable_params(self, directory=None, ):
      """
      Return
      ------
      ([command], stdin, stdout) params for process.Popen
      """

  def __str__(self):
      from io import StringIO
      s = StringIO()
      self.save_to_file(s)
      return s.getvalue()

  def set_option(name, value):
      for i in self._members:
          if name in i:
             self._members[i][name] = value
      else:
          raise KeyError("No option with name {} in any of the members".format(name))

  @staticmethod
  @lazy_value
  def definitions():
      names = (i for i in pkgutil.iter_modules(definitions.__path__) if i.name != 'sections')
      im = importlib.import_module
      modules = ( im('.definitions.' + i.name, __package__) for i in names )
      return { m.task.name.upper(): m.task for m in modules }

  @classmethod
  def is_it_a_task_name(cls, name):
      name = name.upper()
      return name if name in cls.definitions() else False

  @classmethod
  def create_task(cls, arg):
      if isinstance(arg, str):
         name = cls.is_it_a_task_name(arg)
         if name:
            return cls.create(arg), False
         return cls.from_file(arg), True
      if isinstance(arg, io.IOBase):
         return cls.from_file(arg), False
      return arg

  @classmethod
  def create(cls, name):
      return Task(cls.definitions()[name.upper()])

  @classmethod
  def default_task(cls):
      return cls.create('SCF')

  @classmethod
  def from_file(cls, filename):
      definitions = cls.definitions()
      for d in definitions.values():
          try:
             out = d.read_from_file(filename)
             return out
          except Exception as e:
             last = e
      raise last
