"""Density of states (DOS) reader and result."""

import os

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property


class DosResult(TaskResult):
    """Density of states (DOS) results provides access to the computed
    density of states in the DOS output file, using the :py:attr:dos property."""

    @cached_property
    def dos_filename(self):
        """Density of states file name"""
        fname = self.input_parameters.CONTROL.DATASET() + "_DOS.dos"
        if self.directory:
            fname = os.path.join(self.directory, fname)
        return fname

    @cached_property
    def dos(self):
        """The computed density of states."""
        return self._dos_file.parsed_file()

    @cached_property
    def _dos_file(self):
        try:
            return self.files.set_file_type("DOS", "dos")
        except KeyError as exc:
            raise FileNotFoundError("DOS result file is not listed in the SPR-KKR output.") from exc

    @property
    def output_values(self):
        self._dos_file
        return super().output_values


class DosOutputReader(KkrOutputReader):
    """DOS task has no special parser, it just has a special
    result, that allow easy acces to spc output file"""

    result_class = DosResult
    parser_class = DefaultOutputParser

    @classmethod
    def matches_output(cls, output):
        """Recognize legacy DOS output from its ``FILES`` section."""
        lines = output.splitlines()
        for index, line in enumerate(lines):
            if line.strip() != b"FILES:":
                continue
            started = False
            for file_line in lines[index + 1:]:
                file_line = file_line.strip()
                if not file_line:
                    if started:
                        return False
                    continue
                started = True
                name, separator, _ = file_line.partition(b":")
                if separator and name.strip() == b"DOS":
                    return True
            return False
        return False
