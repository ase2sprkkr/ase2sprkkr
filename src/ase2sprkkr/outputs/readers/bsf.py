"""The Bloch spectral functions (BSF) reader and result."""

from ..task_result import TaskResult, KkrOutputReader, OutputFileResultValue
from .default import DefaultOutputParser
from ...common.decorators import cached_property
from ...output_files.output_files import OutputFile


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
        return self.output_values["bsf"]()

    @cached_property
    def output_values(self):
        return {
            "bsf": OutputFileResultValue(
                "Bloch spectral function", "bsf", self.bsf_filename
            )
        }


class BsfOutputReader(KkrOutputReader):
    """BSF task has no special parser, it just has a special
    result, that allow easy acces to BSF output file"""

    result_class = BsfResult
    parser_class = DefaultOutputParser
