"""This module contains classes, used by parsers of the output files"""

import asyncio
import os
import re
import importlib
from . import readers
from ..common.decorators import cached_property, cached_class_property, maybeclassmethod
from ..common.file_utils import FilePath, filename_from_file
from ..potentials.potentials import Potential
from ..input_parameters import input_parameters as input_parameters
from .result_options import create_files_section


class TaskResult:
    """A base class for a result of a runned task (kkrscf executable)"""

    sfn_generated = None
    """Whether this calculation generated rather than reused its SFN file."""

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
            self.output_file = FilePath(output_filename, self.directory)
            self.files.add_file("output", self.output_file)

        input_filename = filename_from_file(input_file, None)
        self.input_file = FilePath(input_filename, self.directory) if input_filename is not None else None

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

    def parse_output_file(self, file, potential=None):
        """Parse a registered output file in the context of this result."""

        from ..output_files.output_files import OutputFile

        if potential is None:
            potential = self.any_potential
        return OutputFile.from_file(
            file.path(),
            try_only=file.file_type,
            potential=potential,
        )

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
            reader = KkrOutputReader.class_for_task(task or "")
            reader = reader(None, None, os.path.dirname(file))
            f.seek(0)
            return reader.read_from_file(f)

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


class KkrOutputReader:
    """Parse SPR-KKR output streams into a TaskResult.

    Execution is handled by :mod:`ase2sprkkr.execution`; this class owns only
    task-specific parsing and result construction.
    """

    def __init__(self, input_parameters, calculator, directory, print_output=False):
        self.input_parameters = input_parameters
        self.calculator = calculator
        self.directory = directory
        self.parser = self.parser_class(print_output)

    @property
    def print_output(self):
        """Whether raw process output should be mirrored to the terminal."""
        return bool(self.parser.print_output)

    def create_result(self, output_file, input_file=None):
        """Create the TaskResult before parsing starts."""
        return self.result_class(
            self.input_parameters,
            self.calculator,
            self.directory,
            output_file=output_file,
            input_file=input_file,
        )

    async def parse_streams(self, stdout, stderr, result):
        """Parse simultaneously growing stdout and stderr streams."""
        return await asyncio.gather(
            self.parser.read_output(stdout, result),
            self.parser.read_error(stderr, result),
        )

    def complete_result(self, result, output, error, return_code):
        """Finalize a result after the execution backend completes."""
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

        parsed_output, parsed_error, parsed_return_code = self.parser.read_from_file(
            output_source, error, [result], return_code
        )
        return self.complete_result(
            result, parsed_output, parsed_error, parsed_return_code
        )

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
