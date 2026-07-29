"""This module contains classes, used by parsers of the output files"""

import os
import re
import importlib
from . import readers
from ..common.decorators import cached_property, cached_class_property, maybeclassmethod
from ..common.process_output_reader import run_coro_sync
from ..common.file_utils import filename_from_file
from ..potentials.potentials import Potential
from ..input_parameters import input_parameters as input_parameters
from .result_options import create_files_section
from pathlib import Path


class TaskResult:
    """A base class for a result of a runned task (kkrscf executable)"""

    def __init__(self, input_parameters, calculator, directory, output_file=None, input_file=None):
        self._input_parameters = input_parameters
        self._calculator = calculator
        self.files = create_files_section(self)
        self._directory = filename_from_file(directory, None, realpath=True)

        output_filename = filename_from_file(output_file, None)
        self.output_file = None
        if output_filename is not None:
            if self._directory is None and os.path.isabs(output_filename):
                self._directory = os.path.dirname(output_filename)
            self.output_file = output_filename
            if os.path.isabs(output_filename):
                output_filename = os.path.relpath(output_filename, self.directory)
                self.output_file = output_filename
            self.files.add_file("output", self.output_file)

        input_filename = filename_from_file(input_file, None)
        if input_filename and os.path.isabs(input_filename):
            input_filename = os.path.relpath(input_filename, self.directory)
        self.input_file = input_filename

    @cached_property
    def directory(self):
        return self._directory or os.getcwd()

    def path_to(self, file):
        """return full path to a given file

        ..doctest::
        >>> t = TaskResult(None, None, "/example")
        >>> t.files.add_file("input", "input.txt")
        <Option input of type String with value = input.txt>
        >>> t.path_to("input")
        '/example/input.txt'
        """
        file = self.files[file]()
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
        if "input" in self.files:
            filename = self.path_to("input")
            if os.path.isfile(filename):
                return filename

    @cached_property
    def potential_filename(self):
        """New (output) potential file name"""
        potfil = self.input_parameters.CONTROL.POTFIL()
        if not potfil:
            potfil = self.files["potential"]() if "potential" in self.files else None
        if not potfil:
            raise ValueError("Please set CONTROL.POTFIL of the input_parameters to read the potential")
        if self.directory:
            potfil = os.path.join(self.directory, potfil)
        return potfil

    @cached_property
    def potential(self):
        """The new (output) potential - that contains the converged charge density etc."""
        return Potential.from_file(self.potential_filename)

    def new_task(self, task):
        out = self._calculator.copy_with_potential(self.potential_filename)
        out.input_parameters = task
        if isinstance(task, str) and task.lower() == "jxc":
            out.input_parameters.set("EMIN", self.last_iteration.energy.EMIN())
        return out

    def complete(self, error, return_code):
        self.error = error
        self.return_code = return_code

    @property
    def atoms(self):
        return self.potential.atoms

    @property
    def output_values(self):
        """Files and other values exposed as the public result summary."""
        return self._ordered_output_files()

    def _ordered_output_files(self):
        """Put task result files before common input/output bookkeeping files."""
        files = self.files.items()
        task_files = {name: value for name, value in files.items() if value.file_type}
        common_files = {name: value for name, value in files.items() if not value.file_type}
        task_files.update(common_files)
        return task_files

    @cached_class_property
    def _match_task_regex(self):
        return re.compile(r" TASK\s+ = ([A-Z]+)\s+\n")

    @classmethod
    def from_file(cls, file):
        def read_output_for(task):
            process = KkrOutputReader.class_for_task(task or "")
            process = process(None, None, os.path.dirname(file))
            f.seek(0)
            return process.read_from_file(f)

        with open(file, "rb") as f:
            raw_out = f.read()
            matches = cls._match_task_regex.search(raw_out.decode("utf8"))
            typ = matches[1] if matches and matches[1] != "NONE" else None
            if typ is None:
                if b"Components of Dzyaloshinski-Moriya vector Dij" in raw_out:
                    typ = "jxc"

            out = read_output_for(typ)
            if typ is None:
                if "DOS" in out.files:
                    with open(file, "rb") as f:
                        out = read_output_for("DOS")
            return out


class KkrProcess:
    def __init__(self, reader, output_parser, coroutine, result, callback=None):
        self.reader = reader
        self.output_parser = output_parser
        self.coroutine = coroutine
        self.result = result
        self.callback = callback

    def run(self):
        result = self.result
        try:
            result, error, return_code = run_coro_sync(self.coroutine)
            self.result.complete(error, return_code)
        finally:
            if self.callback:
                self.callback(self.result)
        return self.result

    def stop_the_process(self):
        self.output_parser.stop_the_process()

    async def run_async(self):
        result, error, return_code = await self.coroutine
        self.result.complete(error, return_code)
        return self.result


class KkrOutputReader:
    """Class, that runs a process and reads its output using the underlying
    output parser (see :class:`ase2sprkkr.common.process_output_reader.ProcessOutputParser`)
    and return the appropriate TaskResult.

    Descendants should define parser_class and result_class property.
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
        self.parser = self.parser_class(print_output, read_callback)

    def _create_result(self, output_file, input_file=None):
        """Create an object that stores results of the KKR output parsing."""
        return self.result_class(
            self.input_parameters, self.calculator, self.directory, output_file=output_file, input_file=input_file
        )

    def _create_process(self, coroutine, result, callback=None):
        """Create an object that can run the desired process (either reading from file or running a process)

        Returns
        -------
        KkrProcess
            An object that can run the process and return the result.
        """
        return KkrProcess(self, self.parser, coroutine, result, callback=callback)

    def create_process(self, cmd, outfile, input_file=None, callback=None, **kwargs):
        """Create an object that takes care of running the command and parsing the results"""
        result = self._create_result(filename_from_file(outfile, None), input_file)
        coroutine = self.parser.run_async(cmd, outfile, self.directory, [result], **kwargs)
        return self._create_process(coroutine, result, callback=callback)

    @maybeclassmethod
    def read_from_file(self, cls, output, error=None, return_code=0, input_file=None, directory=None):
        """Creates an object that takes care of reading and parsing of the output of a sprkkr process"""
        if self is None:
            self = cls(None, None, directory)
        result = self._create_result(output, input_file)
        output_source = output if hasattr(output, "read") else result.path_to("output")

        coroutine = self.parser.read_output_file(output_source, error, [result], return_code)
        return self._create_process(coroutine, result).run()

    @staticmethod
    def class_for_task(task):
        if not task:
            return DefaultOutputReader
        try:
            mod = importlib.import_module(f".{task.lower()}", readers.__name__)
            clsname = task.title() + "OutputReader"
            cls = getattr(mod, clsname)
            if not cls:
                raise Exception(
                    f"Can not determine the class to read the results of task {task}"
                    "No {clsname} class in the module {oo.__name__}.{task}"
                )
        except ModuleNotFoundError:
            cls = DefaultOutputReader

        return cls


from ..outputs.readers.default import DefaultOutputReader  # NOQA
