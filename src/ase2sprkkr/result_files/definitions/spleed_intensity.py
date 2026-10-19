"""Two-column SPLEED I(E) diagram embedded in the native *_SPEC.out log.

SPR-KKR 9.7 SPLOUT writes ENERGY and RINTG using FORMAT(D12.5,2X,D12.5).
The underlying log is not a standalone structured spectrum, so this
ResultFileDefinition extracts its native records before creating a ResultFile.
"""

import os
import re

import numpy as np

from ..result_files import ResultFile
from ..result_files_definitions import (
    ResultFileDefinition,
    ResultFileValueDefinition as V,
)
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV
from ...common.grammar_types import NumpyArray


_DIAGRAM = re.compile(
    r"^\s*(?:i\s*\(e\)\s*-?diagram|rotation diagram|scattering-angle diagram,?)\s*$",
    re.IGNORECASE,
)
# SPLOUT's FORMAT(D12.5,2X,D12.5): exact record widths prevent
# accepting other arbitrary pairs of numbers printed by KKRSPEC.
_NUMBER = re.compile(r"\s*[+-]?\d\.\d{5}[dD][+-]\d{2}\s*\Z")


class SpleedIntensityNotFound(ValueError):
    """The native two-column SPLEED output block is absent."""


def _read_diagram(filename):
    rows = []
    in_diagram = False
    with open(filename, encoding="utf-8", errors="replace") as output:
        for line in output:
            if _DIAGRAM.fullmatch(line.strip()):
                # A new diagram marks a new block, never merge diagrams.
                if rows:
                    break
                in_diagram = True
                continue
            if not in_diagram:
                continue
            record = line.rstrip("\r\n")
            if len(record) != 26 or record[12:14] != "  ":
                continue
            first, second = record[:12], record[14:]
            if _NUMBER.fullmatch(first) and _NUMBER.fullmatch(second):
                rows.append((
                    float(first.replace("D", "E").replace("d", "e")),
                    float(second.replace("D", "E").replace("d", "e")),
                ))
    if not rows:
        raise SpleedIntensityNotFound(
            "No two-column SPLOUT intensity records were found in the SPEC output"
        )
    return np.asarray(rows)


class SpleedIntensityResultFile(ResultFile):
    """Energy and diffracted-beam intensity recorded by SPLOUT in NOUT1."""


class SpleedIntensityDefinition(ResultFileDefinition):
    result_class = SpleedIntensityResultFile
    auto_detect = False

    @property
    def extension(self):
        return "out"

    def read_from_file(self, file, allow_dangerous=False, **kwargs):
        if not isinstance(file, (str, os.PathLike)):
            raise TypeError("SPLEED intensity extraction requires a named SPEC output file")
        filename = os.fspath(file)
        data = _read_diagram(filename)
        result = self.read_from_dict({"RAW_DATA": data}, **kwargs)
        result._filename = filename
        return result


def column(index):
    return slice(None), index


definition = SpleedIntensityDefinition(
    "SPLEED_INTENSITY",
    [
        V("RAW_DATA", NumpyArray(shape=(-1, 2), written_shape=(-1, 2)), name_in_grammar=False),
        NV("ENERGY", "RAW_DATA", column(0), info="Incident electron kinetic energy (eV)"),
        NV("INTENSITY", "RAW_DATA", column(1), info="Diffracted beam intensity (RINTG)"),
    ],
)
