"""This module contains classes used by parsers of SPR-KKR output files."""

import asyncio
import os
import pyparsing as pp
from ..common.warnings import DataValidityError
from ..common.decorators import cached_property, maybeclassmethod
from ..common.file_utils import FilePath, filename_from_file
from ..potentials.potentials import Potential
from ..input_parameters import input_parameters as input_parameters
from .result_options import create_files_section


class TaskResult:
    """A base class for a result of a runned task (kkrscf executable)"""

    produces_potential = False
    """Whether this result produces a new potential usable by later tasks."""

    sfn_generated = None
    """Whether this calculation generated rather than reused its SFN file."""

    def __init__(
        self,
        input_parameters,
        calculator,
        directory,
        output_file=None,
        input_file=None,
        task_name=None,
    ):
        self._input_parameters = input_parameters
        self._detected_input_parameters = None
        self._calculator = calculator
        self._task_name = (
            str(task_name).strip().lower() if task_name is not None else None
        )
        self.files = create_files_section(self)
        self._directory = filename_from_file(directory, None, realpath=True)

        output_filename = filename_from_file(output_file, None)
        self.output_file = None
        if output_filename is not None:
            if self._directory is None and os.path.isabs(output_filename):
                self._directory = os.path.dirname(output_filename)
            self.output_file = FilePath(output_filename, self.directory)
            self.files.add_file("output", self.output_file)

        input_filename = filename_from_file(input_file, None)
        self.input_file = FilePath(input_filename, self.directory) if input_filename is not None else None
        if self.input_file is not None:
            self.files.add_file("input", self.input_file)

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
        return os.fspath(self.files[file]())

    @property
    def task_name(self):
        if self._task_name is not None:
            return self._task_name
        parameters = self.input_parameters
        task = parameters.task_definition if parameters is not None else None
        if task is None:
            raise ValueError("Cannot determine the task of this result from its input")
        return str(task.name).strip().lower()

    @cached_property
    def input_parameters(self):
        if self._input_parameters is not None:
            return self._input_parameters
        source = self.input_parameters_file
        if source is None:
            return None
        if self._detected_input_parameters is not None:
            return self._detected_input_parameters
        # An explicitly registered input is authoritative: parsing errors
        # here must remain visible to the caller.
        return input_parameters.InputParameters.from_file(source)

    def _candidate_input_files(self):
        """Use an unambiguous native input file; never guess between runs."""
        import glob

        directory = self.directory
        if not os.path.isdir(directory):
            return ()
        candidates = sorted(
            {path for suffix in (".inp", ".in")
             for path in glob.glob(os.path.join(directory, "*" + suffix))}
        )
        # A paired input name is authoritative even if the directory contains
        # unrelated calculations.
        if self.output_file is not None:
            stem = os.path.splitext(os.path.basename(os.fspath(self.output_file)))[0]
            for suffix in (".inp", ".in"):
                paired = os.path.join(directory, stem + suffix)
                if paired in candidates:
                    return (paired,)
        return tuple(candidates) if len(candidates) == 1 else ()

    @cached_property
    def input_parameters_file(self):
        if "input" in self.files:
            filename = self.path_to("input")
            if os.path.isfile(filename):
                return filename
            return None
        for filename in self._candidate_input_files():
            try:
                parameters = input_parameters.InputParameters.from_file(filename)
            except (OSError, UnicodeError, ValueError, pp.ParseBaseException, DataValidityError):
                continue
            if self._task_name is not None and parameters.task_name != self._task_name.upper():
                continue
            self._detected_input_parameters = parameters
            self.files.add_file("input", filename)
            return filename
        return None

    @cached_property
    def potential_filename(self):
        """New (output) potential file name"""
        parameters = self.input_parameters
        potfil = parameters.CONTROL.POTFIL() if parameters is not None else None
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

    @property
    def any_potential(self):
        """Return the parsed potential (SCF then try output and then input) """
        try:
            return self.potential
        except (AttributeError, KeyError, OSError, ValueError):
            return None

    def parse_result_file(self, file, potential=None):
        """Parse a registered structured result file in this calculation context."""

        from ..result_files import ResultFile

        if potential is None:
            potential = self.any_potential
        return ResultFile.from_file(
            file.path(),
            try_only=file.file_type,
            potential=potential,
        )

    def new_task(self, task):
        """Create a fresh calculator for *task* using values from this result."""
        from ..input_parameters.task_definition import TaskDefinition
        from ..sprkkr.calculator import SPRKKR

        task = (
            task
            if isinstance(task, TaskDefinition)
            else input_parameters.InputParameters.task_registry[task]
        )
        parameters = (
            input_parameters.InputParameters.input_parameters_registry.create_from_result(
                task, self
            )
        )
        return SPRKKR(
            potential=self.potential_filename,
            input_parameters=parameters,
        )

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

    @classmethod
    def from_file(cls, file):
        """Read a standalone SPR-KKR output using registry-owned detection."""
        from .reader_registry import reader_registry

        filename = os.fspath(file)
        task = reader_registry.detect(filename)
        reader_class = reader_registry.class_for_task_name(task)
        reader = reader_class(
            None,
            None,
            os.path.dirname(filename),
            task_name=task,
        )
        return reader.read_from_file(filename)


class KkrOutputReader:
    """Parse SPR-KKR output streams into a TaskResult.

    Execution is handled by :mod:`ase2sprkkr.execution`; this class owns only
    task-specific parsing and result construction.
    """

    def __init__(
        self,
        input_parameters,
        calculator,
        directory,
        print_output=False,
        task_name=None,
    ):
        self.input_parameters = input_parameters
        self.calculator = calculator
        self.directory = directory
        task = input_parameters.task_definition if input_parameters is not None else None
        self.task_name = (
            task_name
            if task_name is not None
            else task.name
            if task is not None
            else None
        )
        self.parser = self.parser_class(print_output)

    @property
    def print_output(self):
        """Whether raw process output should be mirrored to the terminal."""
        return bool(self.parser.print_output)

    @classmethod
    def matches_output(cls, output):
        """Return whether legacy output without a task marker belongs to this reader."""
        return False

    def create_result(self, output_file, input_file=None):
        """Create the TaskResult before parsing starts."""
        return self.result_class(
            self.input_parameters,
            self.calculator,
            self.directory,
            output_file=output_file,
            input_file=input_file,
            task_name=self.task_name,
        )

    async def parse_streams(self, stdout, stderr, result):
        """Parse simultaneously growing stdout and stderr streams.

        Parser failures are returned to the execution backend rather than
        raised immediately. The backend must first observe the process exit
        status: a truncated output stream is often only a symptom of the
        executable having failed, and that failure is the more useful error.
        """
        return await asyncio.gather(
            self.parser.read_output(stdout, result),
            self.parser.read_error(stderr, result),
            return_exceptions=True,
        )

    def complete_result(self, result, output, error, return_code):
        """Finalize a result after the execution backend completes."""
        if return_code != 0:
            self.parser.result(
                None if isinstance(output, BaseException) else output,
                None if isinstance(error, BaseException) else error,
                return_code,
            )
        if isinstance(output, BaseException):
            raise output
        if isinstance(error, BaseException):
            raise error

        _output, error, return_code = self.parser.result(
            output, error, return_code
        )
        result.complete(error, return_code)
        return result

    @maybeclassmethod
    def read_from_file(self, cls, output, error=None, return_code=0, input_file=None, directory=None):
        """Parse output of an already finished SPR-KKR calculation."""
        if self is None:
            self = cls(None, None, directory)
        result = self.create_result(output, input_file)
        output_source = output if hasattr(output, "read") else result.path_to("output")

        _parsed_output, parsed_error, parsed_return_code = self.parser.read_from_file(
            output_source, error, [result], return_code
        )
        result.complete(parsed_error, parsed_return_code)
        return result
