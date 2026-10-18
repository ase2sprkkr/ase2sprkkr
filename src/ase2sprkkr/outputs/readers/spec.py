"""
Parser for spec.out files containing spectral function calculations.
"""

import re
import numpy as np
from typing import Dict

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property


class SpecResult(TaskResult):
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

    @cached_property
    def spectral_data(self) -> np.ndarray:
        """Extract the main spectral data table."""
        if not hasattr(self, "_spectral_data"):
            self._parse_spectral_data()
        return self._spectral_data

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

    def _parse_spectral_data(self) -> None:
        """Read numeric rows from native KKRSPEC .spc output, not stdout.

        BAND writes eigenvalue rows with a different column count from
        photoemission and SPLEED; preserve all native numeric columns.
        """
        from .arpes import ArpesResult

        filenames = (
            ArpesResult.spc_filenames.func(self)
            if hasattr(self, "_input_parameters") else ()
        )
        if not filenames:
            # Preserve parsing older standalone SPEC text dumps.
            content = self._get_file_content()
            marker = re.search(r"stookes-vector\\s+of the photonfield", content)
            rows = content[marker.end():].splitlines() if marker else ()
        else:
            with open(filenames[0], encoding="utf-8") as file:
                rows = file.readlines()
        data = []
        for line in rows:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            try:
                values = [float(item.replace("D", "E").replace("d", "e"))
                          for item in line.split()]
            except ValueError:
                continue
            if values and (not data or len(values) == len(data[0])):
                data.append(values)
        self._spectral_data = np.asarray(data) if data else np.array([])

    def _get_file_content(self) -> str:
        """Get the content of the output file."""
        if not hasattr(self, "_file_content"):
            if hasattr(self, "output_file") and self.output_file:
                with open(self.output_file, "r") as f:
                    self._file_content = f.read()
            else:
                self._file_content = ""
        return self._file_content


class SpecOutputParser(DefaultOutputParser):
    """Reader for spec.out files."""

    async def read_output(self, stdout, result):
        """Read and parse the output from the spec.out file."""
        await super().read_output(stdout, result)

        # Additional parsing specific to spec.out
        if hasattr(result, "output_lines"):
            result.output_content = "\n".join(result.output_lines)

        return result


class SpecOutputReader(KkrOutputReader):
    """Process class for spectral function calculations."""

    result_class = SpecResult
    parser_class = SpecOutputParser
