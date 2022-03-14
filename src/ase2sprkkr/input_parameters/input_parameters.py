""" Containers for configuration parameters of SPR-KKR task.

These configuration parameters are supplied to SPR-KKR executables as input files.
The containers are also the objects, that take care about executing the SPR-KKR
executables.
"""

import os
import io
import tempfile
import pkgutil
import importlib
from . import definitions
from ..outputs import readers
from ..outputs.readers.default import DefaultOutputReader
from ..common.configuration_containers import RootConfigurationContainer
from ..common.misc import lazy_value, OrderedDict
import shutil
from typing import Union

def resolve_executable_postfix(postfix:Union[str,bool]):
    """" Return the postfix, that is appended after the name of SPR-KKR executable.

    Parameters
    ----------
    postfix
      - If str is given, it is left as is.
      - If True, return the content of SPRKKR_EXECUTABLE_SUFFIX environment variable
      - If False, return ''
    """
    if postfix is False:
        return ''
    if postfix is True:
        return os.getenv('SPRKKR_EXECUTABLE_SUFFIX', '')
    return postfix

class InputParameters(RootConfigurationContainer):
  """ It holds the configuration values for a SPR-KKR task and run the task

  This class is a ConfigurationContainer, thus, it holds the configuration values
  for the task. Moreover, according to its definition, it can run the task -
  execute the executable with proper parameters - and instantiate the result class,
  which then parse the output of the task.
  """

  def __init__(self, definition, inputfile=None, outputfile=False):
      super().__init__(definition)
      self._inputfile = inputfile
      self._outputfile = inputfile if outputfile is False else outputfile

  @property
  def task_name(self):
      """ Return the task name, as defined in the definition of the parameters
      (see InputParametersDefinition class) """
      return self._definition.name

  _default_mpi_runner = None

  @classmethod
  def default_mpi_runner(cls):
      """ Return the executable and its params to run a mpi task.
          The runner is determined by autodetection.

      Return
      ------
      mpi_runner: list or False
          List of strings with executable and its parameters to run a mpi task.
          E.g. [ 'mpirun' ]
          If no suitable runner is found, return False.
      """
      if cls._default_mpi_runner is None:
         for r in [ 'mpirun', 'mpirun.opmpirun', 'mpirun.mpich' ]:
             if shutil.which(r):
                cls._default_mpi_runner = [ r ]
                return cls._default_mpi_runner
         print("No MPI runner found. Disabling MPI!!!")
         cls._default_mpi_runner=False
      return cls._default_mpi_runner

  def mpi_runner(self, mpi):
      """ Return a shell command to execute a mpi task.

      Parameters
      ----------
      mpi_runner: Union[bool|str|list|int]


        - If True is given, return the default mpi-runner
        - If False is given, no mpi-runner is returned.
        - If a string is given, it is interpreted as a list of one item
        - If a list (of strings) is given, the user specified its own runner, use it as is
          as the parameters for subprocess.run.
        - If an integer is given, it is interpreted as the number of
          processes: the default mpi-runner is used, and the parameters
          to specify the number of processes.

      Return
      ------
      mpi_runner: list
        List of strings with the executable and its parameters, e.g.

        ::

            ['mpirun', '-np', '4']
      """
      if mpi is False or not self._definition.mpi:
          return None
      if isinstance(mpi, list):
          return mpi
      if isinstance(mpi, str):
          return [ mpi ]
      runner = self.default_mpi_runner()
      if not runner:
         return None
      if mpi is True:
         return runner
      if isinstance(mpi, int):
         return runner + ['-np', str(mpi)]
      raise Exception("Unknown mpi argument: " + str(mpi))

  def is_mpi(self, mpi=True):
      """ Will be this task (the task described by this input parameters)
      runned using mpi?

      Parameters
      ----------
      mpi: list or str or bool
        Optional parameter. If False is given, return False regardless the type of the
        task, otherwise it is ignored. See InputParameters.mpi_runner for the explanation of the
        parameter.

      Return
      ------
      is_mpi: bool
        Will be this task runned using mpi?
      """
      return mpi and self._definition.mpi

  def run_process(self, calculator, input_file, output_file, print_output=False, executable_postfix=None, mpi=None):
      """
      Run the process that calculate the task specified by this input paameters

      calculator: ase2sprkkr.sprkkr.calculator.SPRKKR
        Calculator, used for running the task. Various configurations are readed from them.

      input_file: file
        File, where the input parameters for the task are stored.

      output_file: file
        File, where the output will be writed. It should be opened for writing (in binary mode)

      print_output: bool or str
        Print the output to the stdout, too? String value 'info' prints only selected infromations
        (depending on the task runned)

      executable_postfix: str or None
        Postfix, appended to the name of the called executable (sometimes, SPRKKR executables are
        compiled so that the resulting executables has postfixies)

      mpi: list or str or bool
        Run the task using mpi? See InputParameters.mpi_runner for the possible values of the parameter.

      Return
      ------
      out: mixed
        The result of ase2sprkkr.common.process_output_reader.ProcessOutputReader.run() method: the
        parsed output of the runned task.

      """
      d = self._definition
      executable = d.executable
      print_output = print_output if print_output is not None else calculator.print_output
      executable_postfix = calculator.executable_postfix if executable_postfix is None else executable_postfix
      executable += resolve_executable_postfix(calculator.executable_postfix)
      self.directory = calculator._directory
      process = self.result_reader(calculator)
      try:
        mpi = self.mpi_runner(mpi)
        if mpi:
             executable = mpi + [ executable + 'MPI', input_file.name ]
             return process.run(executable, output_file, stdin = None, print_output=print_output, directory=self.directory)
        else:
             executable = [ executable ]
             return process.run(executable, output_file, stdin = input_file, print_output=print_output, directory=self.directory)
      except FileNotFoundError as e:
        e.strerror = 'Cannot find SPRKKR executable. Maybe, the SPRKKR_EXECUTABLE_SUFFIX environment variable should be set?\n' + \
                     e.strerror
        raise
      finally:
        input_file.close()

  def result_reader(self, calculator=None):
      """ Return the result readed: the class that parse the output
      of the runned task """
      cls = self._definition.result_reader
      if cls is None:
         task = self.TASK.TASK().lower()
         try:
            mod = importlib.import_module(f'.{task}', readers.__name__)
            clsname = task.title() + 'OutputReader'
            cls = getattr(mod, clsname)
            if not cls:
               raise Exception(f"Can not determine the class to read the results of task {task}"
                                "No {clsname} class in the module {oo.__name__}.{task}")
         except ModuleNotFoundError:
            cls = DefaultOutputReader
      return cls(self, calculator)

  def read_output_from_file(self, filename, directory=None):
      """ Read output of a previously runned task from a file and parse it in a same
      way, as the process would be runned again.
      """
      self.directory = directory or os.path.dirname(filename)
      return self.result_reader().read_from_file(filename)

  def executable_params(self, directory=None, ):
      """
      Return
      ------
      ([executable], stdin, stdout) params for process.Popen
      """

  def __str__(self):
      return f"<Configuration container {self._get_path()}>"

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
  def create_input_parameters(cls, arg):
      """
      Create input_parameters object

      Parameters
      ----------
      arg: str or InputParameters
        If an InputParameters object is given, it is returned as is.
        If a string is given, it is interpreted either as a filename
        (from which the parameters are read) or the task name, for
        which the default parameters are used

      Return
      ------
      input_parameters: InputParameters
      """
      if isinstance(arg, str):
         name = cls.is_it_a_input_parameters_name(arg)
         if name:
            return cls.create(arg)
         return cls.from_file(arg)
      if isinstance(arg, io.IOBase):
         return cls.from_file(arg)
      return arg

  @classmethod
  def create(cls, name):
      """ Create input parameters for the given task name
      Parameters
      ----------
      name: str
        Name of the task (e.g. 'SCF', 'PHAGEN')

      Return
      ------
      input_parameters: InputParameters
        Input parameters with the default values for the given task.
      """
      return InputParameters(cls.definitions()[name.upper()])

  @classmethod
  def default_parameters(cls):
      """ Create default input parameters """
      return cls.create('SCF')

  @classmethod
  def from_file(cls, filename):
      """ Read an input file and create a new InputParameters object from the readed stuff

      Parameters
      ----------
      filename: str or file
        Input file (either its filename, or an open file)
      """
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
      calc = calculator.SPRKKR()
      calc.calculate(input_parameters=self, *args, **kwargs)

#at least, to avoid a circular import
from ..sprkkr import calculator
