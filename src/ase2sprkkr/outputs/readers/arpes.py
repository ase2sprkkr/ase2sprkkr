"""The ARPES reader and result."""

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property
import os
from pathlib import Path


class ArpesResult(TaskResult):
    """ARPES result provides acces to the generated SPC output file
    containing the results of Angle resolved photoelectron electroscopy,
    using the :py:attr:`~spc` property."""

    def to_nexus(self, filename, *, input_file=None, metadata=None):
        """Export ARPES spectrum and recover input provenance when available."""
        if input_file is None:
            input_file = self.input_parameters_file
        else:
            input_file = Path(input_file)
            if not input_file.is_absolute():
                input_file = Path(self.directory) / input_file
        return self.spc.to_nexus(
            filename,
            input_parameters=self._input_parameters,
            input_file=input_file,
            metadata=metadata,
        )

    @cached_property
    def spc_filename(self):
        """Spectroscopy results file name"""
        parameters = self.input_parameters
        if parameters is not None:
            return os.path.join(
                self.directory, parameters.CONTROL.DATASET() + "_ARPES_data.spc"
            )
        matches = list(Path(self.directory).glob("*_ARPES_data.spc"))
        if len(matches) == 1:
            return str(matches[0])
        if not matches:
            raise FileNotFoundError(
                f"No ARPES spectrum found in {self.directory}"
            )
        raise ValueError(
            f"Several ARPES SPC files exist in {self.directory}; "
            "provide a corresponding input file or choose the SPC explicitly"
        )

    @cached_property
    def spc(self):
        """Spectroscopy results."""
        return self._spc_file.parsed_file()

    @property
    def arpes(self):
        """ Alias of ``meth:ArpesResult.spc`` """
        return self.spc

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
