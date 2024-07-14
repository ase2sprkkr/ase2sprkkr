""" Containers for configuration parameters of SPR-KKR task.

These configuration parameters are supplied to SPR-KKR executables as input files.
The containers are also the objects, that take care about executing the SPR-KKR
executables.
"""
from __future__ import annotations

import os
import warnings
import io
import pkgutil
import importlib
from . import definitions
from ..outputs import readers
from ..outputs.readers.default import DefaultProcess
from ..sprkkr.configuration import ConfigurationFile, ConfigurationSection
from ..common.decorators import cached_class_property
import shutil
from typing import Union
from .. import config


class InputSection(ConfigurationSection):
    """ Input parameters sections has nothing special, yet. """


class InputParameters(ConfigurationFile):
  """ It holds the configuration values for a SPR-KKR task and run the task

  This class is a ConfigurationContainer, thus, it holds the configuration values
  for the task. Moreover, according to its definition, it can run the task -
  execute the executable with proper parameters - and instantiate the result class,
  which then parse the output of the task.
  """

  def resolve_executable_suffix(self, postfix:Union[str,bool]):
      """" Return the postfix, that is appended after the name of SPR-KKR executable.

      Parameters
      ----------
      postfix
        - If str is given, it is left as is.
        - If True, return the default value:
                    config.sprkkr_executable_suffix
                   (content of SPRKKR_EXECUTABLE_SUFFIX env variable
                    or a user_defined_value)
        - If False, return ''
      """
      if postfix is False:
          return ''
      if postfix is True:
          return config.sprkkr_executable_suffix
      return postfix

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
  def default_mpi_runner(cls, auto:bool=False):
      """ Return the executable and its params to run a mpi task.
          The runner is determined by autodetection.

      Parameters
      ----------
      auto:bool
          If true,
             - do not warn if no mpi is found
             - check the number of available cpus, do not return a runner, if there is only one

      Return
      ------
      mpi_runner: list or False
          List of strings with executable and its parameters to run a mpi task.
          E.g. [ 'mpirun' ]
          If no suitable runner is found, return False.
      """
      if auto and hasattr(os, 'sched_getaffinity') and len(os.sched_getaffinity(0))==1:
          return None

      if config.mpi_runner:
          return config.mpi_runner

      if cls._default_mpi_runner is None:
         for r in [ 'mpirun', 'mpirun.opmpirun', 'mpirun.mpich' ]:
             if shutil.which(r):
                cls._default_mpi_runner = [ r ]
                return cls._default_mpi_runner
         if not auto:
             warnings.warn("No MPI runner found. Disabling MPI!!!")
         cls._default_mpi_runner=False
      return cls._default_mpi_runner

  def mpi_runner(self, mpi):
      """ Return a shell command to execute a mpi task.

      Parameters
      ----------
      mpi_runner: Union[bool,str,list,int]


        - If True is given, return the default mpi-runner
        - If False is given, no mpi-runner is returned.
        - If 'auto' is given, it is the same as True, however
             * no warning is given if no mpi is found
             * MPI is not used, if only one CPU is available
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
          if mpi == 'auto':
              return self.default_mpi_runner(auto=True) or None
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

  def run_process(self, calculator, input_file, output_file, directory='.', print_output=None, executable_suffix=None, mpi=None, gdb=False):
      """
      Run the process that calculate the task specified by this input paameters

      calculator: ase2sprkkr.sprkkr.calculator.SPRKKR
        Calculator, used for running the task. Various configurations are readed from them.

      directory: str
        Where to run the calculation

      input_file: file
        File, where the input parameters for the task are stored.

      output_file: file
        File, where the output will be writed. It should be opened for writing (in binary mode)

      print_output: bool or str
        Print the output to the stdout, too? String value 'info' prints only selected infromations
        (depending on the task runned)

      executable_suffix: str or None
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
      print_output = calculator.value_or_default('print_output', print_output)
      executable_suffix = calculator.value_or_default('executable_suffix', executable_suffix)
      mpi = calculator.value_or_default('mpi', mpi)

      d = self._definition
      executable = d.executable
      executable += self.resolve_executable_suffix(executable_suffix)
      process = self.result_reader(calculator)
      try:
        mpi = self.mpi_runner(mpi)
        if mpi:
             executable = mpi + [ executable + 'MPI', input_file.name ]
             stdin = None
        else:
             executable = [ executable ]
             stdin = input_file
        if gdb:
             stdin = None
             if mpi:
                print(f"run {input_file.name}")
                executable = executable[:-1]
             else:
                print(f"run <{input_file.name}")
             executable = [ 'gdb' ] + executable

        return process.run(executable, output_file, stdin = stdin, print_output=print_output, directory=directory, input_file=input_file.name)
      except FileNotFoundError as e:
        e.strerror = 'Cannot find SPRKKR executable. Maybe, the SPRKKR_EXECUTABLE_SUFFIX environment variable or ase2sprkkr.config.sprkkr_executable_suffix variable should be set?\n' + \
                     e.strerror
        raise
      finally:
        input_file.close()

  def result_reader(self, calculator=None, directory=None):
      """ Return the result readed: the class that parse the output
      of the runned task

      calculator
        Calculator, which will be attached to the resulting class
        for ruther processing the results

      directory
        Directory, to which will be related the relative paths
        in the result.
        If none, get the directory from the calculator, or the current
        directory
      """
      cls = self._definition.result_reader

      if not directory and calculator.directory:
         directory = calculator.directory

      if cls is None:
         task = self.TASK.TASK().lower()
         try:
            mod = importlib.import_module(f'.{task}', readers.__name__)
            clsname = task.title() + 'Process'
            cls = getattr(mod, clsname)
            if not cls:
               raise Exception(f"Can not determine the class to read the results of task {task}"
                                "No {clsname} class in the module {oo.__name__}.{task}")
         except ModuleNotFoundError:
            cls = DefaultProcess
      return cls(self, calculator, directory)

  def read_output_from_file(self, filename, directory=None):
      """ Read output of a previously runned task from a file and parse it in a same
      way, as the process would be runned again.

      filename
        Filename, from which the result will be read

      directory
        A directory, to which are related the paths in the output
        Default None means the directory, where the file is
      """
      directory = directory or os.path.dirname(filename)
      return self.result_reader(directory=directory).read_from_file(filename)

  def executable_params(self, directory=None, ):
      """
      Return
      ------
      ([executable], stdin, stdout) params for process.Popen
      """

  def __str__(self):
      return f"<Configuration container {self._get_path()}>"

  def set_option(self, name, value):
      for i in self._members:
          if name in i:
             self._members[i][name] = value
      else:
          raise KeyError("No option with name {} in any of the members".format(name))

  @cached_class_property
  def definitions():
      # user = os.path.join(platformdirs.user_config_dir('ase2sprkkr', 'ase2sprkkr'), 'input_parameters')
      names = (i for i in pkgutil.iter_modules(definitions.__path__))
      im = importlib.import_module
      modules = ( im('.definitions.' + i.name, __package__) for i in names )
      return { m.__name__.rsplit('.',1)[1].upper(): m.input_parameters for m in modules if hasattr(m, 'input_parameters') }

  @classmethod
  def is_it_a_input_parameters_name(cls, name):
      name = name.upper()
      return name if name in cls.definitions else False

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
      return InputParameters(cls.task_definition(name))

  @classmethod
  def task_definition(cls, task_name):
      return cls.definitions[task_name.upper()]

  @classmethod
  def default_parameters(cls):
      """ Create default input parameters """
      return cls.create('SCF')

  @classmethod
  def from_file(cls, filename, allow_dangerous:bool=False):
      """ Read an input file and create a new InputParameters object from the readed stuff

      Parameters
      ----------
      filename: str or file
        Input file (either its filename, or an open file)

      allow_dangerous
        Allow to read dangerous values of options: i.e. the values that do not fullfil the
        required type of the given option or its other requirements.
      """
      definitions = cls.definitions
      for d in definitions.values():
          try:
             out = d.read_from_file(filename, allow_dangerous=allow_dangerous)
             return out
          except Exception as e:
             last = e
      raise last

  def calculate(self, *args, **kwargs):
      """ Create a calculator and run the input_parameters. See SPRKKR.calculate for the arguments """
      calc = calculator.SPRKKR()
      calc.calculate(input_parameters=self, *args, **kwargs)

  def __repr__(self):
      d = self._definition
      out = d.configuration_type_name
      out = out + ' for task ' + d.name.upper()
      return out

  def change_task(self, task):
      """ Change the task to the given task. Retain the value of the options,
      that are present in the new task.
      """
      vals = self.to_dict()
      self._definition = self.task_definition(task)
      self._init_members_from_the_definition()
      self.set(vals, unknown = 'ignore', error='ignore')

  def save_to_file(self, file, atoms=None, *, validate='save'):
      if self._definition.save_hook:
          self._definition.save_hook(
             getattr(file, "name", None),
             atoms
          )
      super().save_to_file(file, atoms, validate=validate)


# at least, to avoid a circular import
from ..sprkkr import calculator   # NOQA: E402
