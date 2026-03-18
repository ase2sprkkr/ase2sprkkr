""" This module contains classes, used by parsers of the output files """
import os
import re
import importlib
from . import readers
from ..common.decorators import cached_property, cached_class_property
from ..common.process_output_reader import run_coro_sync
from ..potentials.potentials import Potential
from ..input_parameters import input_parameters as input_parameters
from ..output_files.output_files import OutputFile
from pathlib import Path
import os
import platform
import subprocess
import shutil

class ResultValue:

  def __init__(self, name, value):
      self.name = name
      self._value = value

  def __call__(self):
      return self._value

  def value_label(self):
      return self._value

  def actions(self):
      return []

class OutputFileResultValue:

  def __init__(self, name, type, filename):
      self.name = name
      self.type = type
      self.filename = filename

  def __call__(self):
      return OutputFile.from_file(self.filename, try_only=self.type)

  def definition(self):
      return OutputFile.definitions[self.type]

  def value_label(self):
      return f'<{self.type} file>'

  def actions(self):
      out = [ 'open', 'save' ]
      if self.definition().result_class.can_be_plotted():
          out.append('plot')
      return out

  def save(self, parent=None):
      try:
        from PyQT6.QtWidgets import QFileDialog
      except ImportError:
        raise Exception("PyQt6 is required for saving the output file")

      file_path, _ = QFileDialog.getSaveFileName(
          parent,
          "Save File",
          "",
          f" (*.{self.definition().extension});;All Files (*)"
      )
      shutil.copy(self.filename(), file_path);

  def plot(self, parent=None):
      self().plot()

  def open(self, parent=None):
      if platform.system() == "Windows":
          os.startfile(self.filename)
      elif platform.system() == "Darwin":  # macOS
          subprocess.run(["open", self.filename])
      else:  # Linux
          subprocess.run(["xdg-open", self.filename])


class TaskResult:
  """ A base class for a result of a runned task (kkrscf executable) """
  def __init__(self, input_parameters, calculator, directory,
                     output_file=None, input_file=None):

      def file_name(f):
          if isinstance(f, str):
              return f
          if hasattr(f, 'name'):
              return f.name
          return f

      self._input_parameters = input_parameters
      self._calculator = calculator
      self.output_file = output_file
      self.files={}
      if output_file:
          self.files['output'] = output_file
      self.directory = directory or os.path.dirname(file_name(self.files.get('output')) or '') or os.getcwd()
      self.directory = os.path.realpath(self.directory)
      self.input_file = input_file

  def path_to(self, file):
      """ return full path to a given file

      ..doctest::
      >>> t = TaskResult(None, None, '/example')
      >>> t.files['input'] = 'input.txt'
      >>> t.path_to('input')
      '/example/input.txt'
      """
      file = self.files[file]
      if Path(file).is_absolute():
          return file
      return os.path.join(self.directory, file)

  @property
  def task_name(self):
      return self.input_parameters.task_name

  @cached_property
  def input_parameters(self):
      if self._input_parameters:
          return self._input_parameters
      if self.input_parameters_file:
          return input_parameters.InputParameters.from_file(self.input_parameters_file)

  @cached_property
  def input_parameters_file(self):
      if 'input' in self.files and os.path.isfile(self.files['input']):
          return self.files['input']

  @cached_property
  def potential_filename(self):
      """ New (output) potential file name """
      potfil = self.input_parameters.CONTROL.POTFIL()
      if not potfil:
          potfil = self.files.get('potential', None)
      if not potfil:
         raise ValueError("Please set CONTROL.POTFIL of the input_parameters to read the potential")
      if self.directory:
         potfil = os.path.join(self.directory, potfil)
      return potfil

  @cached_property
  def potential(self):
      """ The new (output) potential - that contains the converged charge density etc. """
      return Potential.from_file(self.potential_filename)

  def new_task(self, task):
      out = self._calculator.copy_with_potential(self.potential_filename)
      out.input_parameters = task
      if isinstance(task, str) and task.lower() == 'jxc':
          out.input_parameters.set('EMIN', self.last_iteration.energy.EMIN())
      return out

  def complete(self, error, return_code):
      self.error = error
      self.return_code = return_code

  @property
  def atoms(self):
      return self.potential.atoms

  @cached_class_property
  def _match_task_regex(self):
      return re.compile(r" TASK\s+ = ([A-Z]+)\s+\n")

  @classmethod
  def from_file(cls, file):
      def read_output_for(task):
          process = KkrProcessRunner.class_for_task(task)
          process = process(None, None, os.path.dirname(file))
          f.seek(0)
          return process.read_from_file(f)

      with open(file, "rb") as f:
          raw_out = f.read()
          matches = cls._match_task_regex.search(raw_out.decode('utf8'))
          out = read_output_for(matches[1])
          if matches[1] == 'NONE':
              if 'DOS' in out.files:
                    with open(file, "rb") as f:
                        out = read_output_for('DOS')
          return out

class KkrProcess:

  def __init__(self, runner, output_reader, coroutine, callback=None):
      self.runner = runner
      self.output_reader = output_reader
      self.coroutine = coroutine
      self.callback = callback

  def run(self):
      result = None
      try:
        result, error, return_code = run_coro_sync(self.coroutine)
        result.complete(error, return_code)
      finally:
        if self.callback:
            self.callback(result)
      return result

  def stop_the_process(self):
       self.output_reader.stop_the_process()

  async def run_async(self):
      await self.coroutine
      result.complete(error, return_code)
      return result

class KkrProcessRunner:
  """ Class, that run a process and read its output using underlined
  process reader (see :class:`ase2sprkkr.common.process_output_reader.ProcessOutputReader`)
  and return the appropriate TaskResult.

  Descendants should define reader_class and result_class property.
  """

  def __init__(self, input_parameters, calculator, directory, print_output=False, read_callback=None):
      self.input_parameters = input_parameters
      """ Input parameters, that command to read the output (thus probably the ones, that
      run the process that produced the output. It is used e.g. for determining the potential file,
      which belongs to the output.
      """
      self.calculator = calculator
      self.directory = directory
      """ Calculator, that can be used for further processing of the results. """
      self.reader = self.reader_class(print_output, read_callback)

  def _create_result(self, output_file, input_file=None):
      """ Create an object that stores results of the KKR output parsing. """
      return self.result_class(self.input_parameters, self.calculator, self.directory,
                               output_file = output_file,
                               input_file = input_file
                               )

  def _create_process(self, coroutine, result,callback=None):
      """ Create an object that can run the desired process (either reading from file or running a process) """
      return KkrProcess(self, self.reader, coroutine, callback=callback)

  def create_process(self, cmd, outfile, input_file=None, callback=None, **kwargs):
      """ Create an object that takes care of running the command and parsing the results """
      result = self._create_result(getattr(outfile, "name", None), input_file)
      coroutine = self.reader.run_async(cmd, outfile, self.directory, [result], **kwargs)
      return self._create_process(coroutine, result, callback=callback)

  def read_from_file(self, output, error=None, return_code=0, input_file=None):
      """ Creates an object that takes care of reading and parsing of the output of a sprkkr process """
      result = self._create_result(output, input_file)
      coroutine = self.reader.read_output_file(output, error, [result], return_code)
      return self._create_process(coroutine, result).run()

  @staticmethod
  def class_for_task(task):
       try:
          mod = importlib.import_module(f'.{task.lower()}', readers.__name__)
          clsname = task.title() + 'ProcessRunner'
          cls = getattr(mod, clsname)
          if not cls:
             raise Exception(f"Can not determine the class to read the results of task {task}"
                              "No {clsname} class in the module {oo.__name__}.{task}")
       except ModuleNotFoundError:
          cls = DefaultProcessRunner

       return cls


from ..outputs.readers.default import DefaultProcessRunner   # NOQA
