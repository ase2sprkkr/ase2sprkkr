import os
import io
import tempfile
import pkgutil
import importlib
from . import definitions
from ..common.conf_containers import RootConfContainer
from ..common.misc import lazy_value, OrderedDict

def resolve_executable_postfix(postfix):
    if postfix is False:
        return ''
    if postfix is True:
        return os.getenv('SPRKKR_EXECUTABLE_SUFFIX', '')
    return postfix

class InputParameters(RootConfContainer):
  def __init__(self, definition, inputfile=None, outputfile=False):
      super().__init__(definition)
      self._inputfile = inputfile
      self._outputfile = inputfile if outputfile is False else outputfile

  @property
  def task_name(self):
      return self._definition.name

  _default_mpi_runner = None

  @classmethod
  def default_mpi_runner(cls):
      """ Return the executable to run mpi obtained by an autodetection """
      if cls._default_mpi_runner is None:
         cls._default_mpi_runner = [ 'mpirun' ] if shutil.which('mpirun') else []
      return cls._default_mpi_runner

  @classmethod
  def mpi_runner(cls, mpi):
        if mpi is True:
           return cls.default_mpi_runner
        return [ mpi ] if isinstance(mpi, str) else mpi

  def run_process(self, calculator, task_file, output_file, print_output=False, executable_postfix=None, mpi=None):
      d = self._definition
      executable = d.executable
      print_output = print_output if print_output is not None else calculator.print_output
      executable_postfix = calculator.executable_postfix if executable_postfix is None else executable_postfix
      executable += resolve_executable_postfix(calculator.executable_postfix)
      if d.mpi and mpi:
           executable = self.mpi_runner(mpi) + [ executable + 'MPI' ]
      else:
           executable = [ executable ]
      self.directory = calculator._directory
      process = self.result_reader(calculator)
      try:
        return process.run(executable, output_file, stdin = task_file, print_output=print_output, directory=self.directory)
      except FileNotFoundError as e:
        e.strerror = 'Cannot find SPRKKR executable. Maybe, the SPRKKR_EXECUTABLE_SUFFIX environment variable should be set?\n' + \
                     e.strerror
        raise

  def result_reader(self, calculator=None):
      return self._definition.result_reader(self, calculator)

  def read_output_from_file(self, filename, directory=None):
      self.directory = directory or os.path.dirname(filename)
      return self.result_reader().read_from_file(filename)

  def executable_params(self, directory=None, ):
      """
      Return
      ------
      ([executable], stdin, stdout) params for process.Popen
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
      return { m.input_parameters.name.upper(): m.input_parameters for m in modules }

  @classmethod
  def is_it_a_input_parameters_name(cls, name):
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
      return InputParameters(cls.definitions()[name.upper()])

  @classmethod
  def default_parameters(cls):
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

  def calculate(self, *args, **kwargs):
      """ Create a calculator and run the input_parameters. See SPRKKR.calculate for the arguments """
      calculator = SPRKKR()
      calculator.calculate(task = task, *args, **kwargs)