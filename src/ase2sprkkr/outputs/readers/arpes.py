"""The ARPES reader and result."""

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property
import os


class ArpesResult(TaskResult):
    """ARPES result provides acces to the generated SPC output file
    containing the results of Angle resolved photoelectron electroscopy,
    using the :py:attr:`~spc` property."""

    @cached_property
    def spc_filename(self):
        """Spectroscopy results file name"""
        fname = self.input_parameters.CONTROL.DATASET() + "_ARPES_data.spc"
        if self.directory:
            fname = os.path.join(self.directory, fname)
        return fname

    @cached_property
    def spc(self):
        """Spectroscopy results."""
        return self._spc_file.parsed_file()

    @property
    def arpes(self):
        """ Alias of ``meth:ArpesResult.spc`` """
        return self.spc()

    @property
    def output_values(self):
        self._spc_file
        return super().output_values

    @cached_property
    def _spc_file(self):
        if "SPC" not in self.files:
            return self.files.add_file("SPC", self.spc_filename, "spc")
        return self.files.set_file_type("SPC", "spc")


class ArpesOutputReader(KkrOutputReader):
    """ARPES task has no special parser, it just has a special
    result, that allow easy acces to spc output file"""

    result_class = ArpesResult
    parser_class = DefaultOutputParser
