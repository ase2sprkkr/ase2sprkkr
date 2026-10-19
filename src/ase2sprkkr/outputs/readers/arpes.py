"""The ARPES reader and result."""

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property
import os
import warnings


class ArpesResult(TaskResult):
    """ARPES result provides acces to the generated SPC output file
    containing the results of Angle resolved photoelectron electroscopy,
    using the :py:attr:`~spc` property."""

    @cached_property
    def spc_filenames(self):
        """Native SPC outputs, declared by the task or discovered alongside stdout."""
        parameters = self.input_parameters
        if parameters is not None:
            return tuple(
                os.fspath(path) for path in parameters.output_files(self.directory)
                if os.fspath(path).endswith(".spc")
                and (not getattr(path, "optional", False) or os.path.isfile(path))
            )
        if self.output_file is None:
            return ()
        output = os.path.basename(os.fspath(self.output_file))
        if not output.endswith("_SPEC.out"):
            return ()
        dataset = output[:-len("_SPEC.out")]
        candidates = (
            os.path.join(self.directory, dataset + suffix)
            for suffix in ("_data.spc", "_data_pol1.spc", "_data_pol2.spc")
        )
        found = tuple(path for path in candidates if os.path.isfile(path))
        return found

    @cached_property
    def spc_filename(self):
        """First SPC output, retained as a convenience for older callers."""
        return next(iter(self.spc_filenames), None)

    @cached_property
    def spectra(self):
        """Tuple of parsed spectral files, including all photon polarizations.

        A single-polarization calculation still returns a one-element tuple.
        """
        return tuple(file.parsed_file() for file in self._register_spc_files())

    @property
    def spc(self):
        """Deprecated single-spectrum accessor; prefer :attr:`spectra`."""
        warnings.warn(
            "ArpesResult.spc is deprecated because a calculation can produce "
            "multiple polarization spectra. Use result.spectra (always a tuple) "
            "and select a spectrum explicitly, e.g. result.spectra[0].",
            DeprecationWarning,
            stacklevel=2,
        )
        if not self.spc_filenames:
            raise FileNotFoundError("No SPC output file could be identified")
        return self.spectra[0]

    @property
    def arpes(self):
        """Legacy alias for the first spectral file."""
        return self.spc

    @property
    def output_values(self):
        self._register_spc_files()
        return super().output_values

    def _register_spc_files(self):
        filenames = self.spc_filenames
        if not filenames:
            return ()
        entries = []
        task = self._task_name or (self.input_parameters.task_name if self.input_parameters is not None else None)
        file_type = {"band": "band_spc", "spleed": "spleed_spc"}.get(
            str(task).lower(), "spc"
        )
        for index, filename in enumerate(filenames):
            key = "SPC" if index == 0 else f"SPC_{index + 1}"
            if key not in self.files:
                entries.append(self.files.add_file(key, filename, file_type))
            else:
                entries.append(self.files.set_file_type(key, file_type))
        return tuple(entries)

    @cached_property
    def _spc_file(self):
        files = self._register_spc_files()
        if not files:
            raise FileNotFoundError("No SPC output file could be identified")
        return files[0]


class ArpesOutputReader(KkrOutputReader):
    """ARPES task has no special parser, it just has a special
    result, that allow easy acces to spc output file"""

    result_class = ArpesResult
    parser_class = DefaultOutputParser
