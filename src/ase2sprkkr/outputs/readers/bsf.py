"""The Bloch spectral functions (BSF) reader and result."""

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property


class BsfResult(TaskResult):
    """BSF result provides access to the computed Bloch spectral functions in the
    BSF output file using the :py:attr:`~bsf` property."""

    @cached_property
    def bsf_filename(self):
        """Bloch spectral functions file name"""
        return self.path_to("Bloch-SF")

    @cached_property
    def bsf(self):
        """The computed Bloch spectral functions."""
        return self._bsf_file.parsed_file()

    @cached_property
    def _bsf_file(self):
        try:
            return self.files.set_file_type("Bloch-SF", "bsf")
        except KeyError as exc:
            raise FileNotFoundError("BSF result file is not listed in the SPR-KKR output.") from exc

    @property
    def output_values(self):
        self._bsf_file
        return super().output_values


class BsfOutputReader(KkrOutputReader):
    """BSF task has no special parser, it just has a special
    result, that allow easy acces to BSF output file"""

    result_class = BsfResult
    parser_class = DefaultOutputParser
