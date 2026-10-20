"""
Parser for spec.out files containing spectral function calculations.
"""

import os
import re
import numpy as np
from typing import Dict

from ..task_result import KkrOutputReader
from .arpes import ArpesResult
from .default import DefaultOutputParser
from ...common.decorators import cached_property


class SpleedIntensity:
    """Numerical I(E) data embedded in a SPLEED calculation's SPEC.out log.

    This is not a separate file. It preserves the familiar callable column
    accessors without routing embedded output through ResultFile.
    """

    def __init__(self, rows):
        self._data = np.asarray(rows, dtype=float)

    def RAW_DATA(self):
        return self._data

    def ENERGY(self):
        return self._data[:, 0]

    def INTENSITY(self):
        return self._data[:, 1]


class SpecResult(ArpesResult):
    """
    Class representing the results of a spectral function calculation.
    Provides access to the parsed data from spec.out files.
    """

    @cached_property
    def potential_barrier(self) -> Dict[str, float]:
        """Extract potential barrier parameters from the output."""
        if not hasattr(self, "_potential_barrier"):
            self._parse_potential_barrier()
        return self._potential_barrier

    @cached_property
    def lattice_constants(self) -> Dict[str, float]:
        """Extract lattice constants from the output."""
        if not hasattr(self, "_lattice_constants"):
            self._parse_lattice_info()
        return self._lattice_constants

    @cached_property
    def basis_vectors(self) -> Dict[str, np.ndarray]:
        """Extract basis vectors from the output."""
        if not hasattr(self, "_basis_vectors"):
            self._parse_lattice_info()
        return self._basis_vectors

    @property
    def spleed_intensity(self):
        """Current two-column SPLEED diagram from the main text output, if any.

        Do not cache: callers may inspect a result while SPEC.out is growing.
        """
        if self.task_name != "spleed" or self.output_file is None:
            return None
        filename = self.path_to("output")
        if not os.path.isfile(filename):
            return None
        with open(filename, encoding="utf-8", errors="replace") as stream:
            return SpecOutputParser.parse_spleed_intensity(stream)

    @property
    def spc_filenames(self):
        """SPLEED's header-only SPC is not a separate spectrum."""
        if self.task_name == "spleed":
            parameters = self.input_parameters
            if parameters is not None and str(parameters.SPEC_EL.TYP()) == "0":
                return ()
            if self.spleed_intensity is not None:
                return ()
        return super().spc_filenames

    @property
    def spectral_data(self) -> np.ndarray:
        """Convenience array from the unique spectrum or legacy text log."""
        if self.output_file is None:
            return SpecOutputParser.parse_inline_spectral_data(
                self._get_file_content()
            )
        if self.task_name == "spleed":
            spectrum = self.spleed_intensity
            if spectrum is not None:
                return spectrum.RAW_DATA()
        spectra = self.spectra
        if len(spectra) > 1:
            raise ValueError("Multiple spectral files; select one from result.spectra")
        if spectra:
            return spectra[0].RAW_DATA()
        if os.path.isfile(self.path_to("output")):
            return SpecOutputParser.parse_inline_spectral_data(
                self._get_file_content()
            )
        return np.array([])

    def _parse_potential_barrier(self) -> None:
        """Parse the potential barrier section of the output."""
        self._potential_barrier = {}
        content = self._get_file_content()

        # Parse potential barrier parameters
        ibar_match = re.search(r"#\s*ibar:\s*(\d+)", content)
        if ibar_match:
            self._potential_barrier["ibar"] = int(ibar_match.group(1))

        epsx_match = re.search(r"#\s*epsx:\s*([\d.E+-]+)", content)
        if epsx_match:
            self._potential_barrier["epsx"] = float(epsx_match.group(1))

        # Parse zparup and zpardn (these appear to be arrays)
        zpar_match = re.search(r"#\s*zparup:\s*([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)", content)
        if zpar_match:
            self._potential_barrier["zparup"] = [float(x) for x in zpar_match.groups()]

        zpardn_match = re.search(r"#\s*zpardn:\s*([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)", content)
        if zpardn_match:
            self._potential_barrier["zpardn"] = [float(x) for x in zpardn_match.groups()]

        # Parse bparp
        bparp_match = re.search(r"#\s*bparp:\s*([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)", content)
        if bparp_match:
            self._potential_barrier["bparp"] = [float(x) for x in bparp_match.groups()]

    def _parse_lattice_info(self) -> None:
        """Parse lattice constants and basis vectors from the output."""
        content = self._get_file_content()
        self._lattice_constants = {}
        self._basis_vectors = {}

        # Parse lattice constant
        lat_match = re.search(r"lattice constant\s*:\s*([\d.]+)", content)
        if lat_match:
            self._lattice_constants["a"] = float(lat_match.group(1))

        self._basis_vectors = {"real": np.zeros((2, 2)), "reciprocal": np.zeros((2, 2))}
        # Parse real basis vectors
        matches = re.finditer(r"(real|reciprocal) basis\s+(\d+)\s*:\s*([\d.-]+)\s+([\d.-]+)", content)
        for match in matches:
            self._basis_vectors[match.group(1)][int(match.group(2)) - 1] = [match.group(3), match.group(4)]

    def _get_file_content(self) -> str:
        """Reuse captured SPEC text; fall back to disk for unparsed results."""
        if not hasattr(self, "_file_content"):
            if hasattr(self, "output_lines"):
                self._file_content = "\n".join(self.output_lines)
            elif self.output_file is not None:
                with open(self.output_file, "r") as stream:
                    self._file_content = stream.read()
            else:
                self._file_content = ""
        return self._file_content


class SpecOutputParser(DefaultOutputParser):
    """Reader for the main text output produced by KKRSPEC."""

    _DIAGRAM = re.compile(
        r"^\s*(?:i\s*\(e\)\s*-?diagram|rotation diagram|"
        r"scattering-angle diagram,?)\s*$",
        re.IGNORECASE,
    )
    _D_NUMBER = re.compile(r"\s*[+-]?\d\.\d{5}[dD][+-]\d{2}\s*\Z")

    @classmethod
    def parse_spleed_intensity(cls, lines):
        """Extract the first SPLOUT D12.5,2X,D12.5 table from SPEC.out.

        The diagram is part of NOUT1, not another result file. Stop at the
        first non-record after data begins; do not collect later diagnostics.
        Return None when only a seven-column diagram or no table is present.
        """
        in_diagram = False
        rows = []
        for line in lines:
            if cls._DIAGRAM.fullmatch(line.strip()):
                if rows:
                    break
                in_diagram = True
                continue
            if not in_diagram:
                continue
            record = line.rstrip("\r\n")
            if len(record) == 26 and record[12:14] == "  ":
                first, second = record[:12], record[14:]
                if cls._D_NUMBER.fullmatch(first) and cls._D_NUMBER.fullmatch(second):
                    rows.append((
                        float(first.replace("D", "E").replace("d", "e")),
                        float(second.replace("D", "E").replace("d", "e")),
                    ))
                    continue
            if record.strip():
                # Blank lines follow SPLOUT's diagram header; a different
                # nonempty line marks the end of the current diagram.
                in_diagram = False
                if rows:
                    break
        return SpleedIntensity(rows) if rows else None

    @staticmethod
    def parse_inline_spectral_data(content):
        """Read legacy photoemission diagnostic rows from a text-only log."""
        marker = re.search(r"stookes-vector\s+of the photonfield", content)
        if marker is None:
            return np.array([])
        rows = []
        for line in content[marker.end():].splitlines():
            try:
                values = [
                    float(item.replace("D", "E").replace("d", "e"))
                    for item in line.split()
                ]
            except ValueError:
                if rows and line.strip():
                    break
                continue
            if values:
                if rows and len(values) != len(rows[0]):
                    break
                rows.append(values)
        return np.asarray(rows) if rows else np.array([])

class SpecOutputReader(KkrOutputReader):
    """Process class for spectral function calculations."""

    result_class = SpecResult
    parser_class = SpecOutputParser
