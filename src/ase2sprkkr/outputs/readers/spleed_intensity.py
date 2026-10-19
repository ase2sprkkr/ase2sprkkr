"""Read the two-column SPLEED intensity records embedded in KKRSPEC output.

In SPR-KKR 9.7, SPLOUT writes ENERGY and RINTG to NOUT1 using FORMAT
(D12.5,2X,D12.5). This is not the seven-column SPLEED SPC format.
"""

import re
import numpy as np

_DIAGRAM = re.compile(
    r"\b(?:i\s*\(e\)\s*-?diagram|rotation diagram|scattering-angle diagram)",
    re.IGNORECASE,
)
_D_NUMBER = re.compile(r"[+-]?(?:\d+\.\d*|\.\d+)[dD][+-]?\d{2,3}\Z")


class SpleedIntensity:
    """SPLOUT two-column data: incident energy (eV), diffracted intensity."""

    def __init__(self, data):
        values = np.asarray(data, dtype=float)
        if values.ndim != 2 or values.shape[1] != 2:
            raise ValueError("SPLEED intensity data require exactly two columns")
        self._data = values

    def RAW_DATA(self):
        return self._data

    def ENERGY(self):
        return self._data[:, 0]

    def INTENSITY(self):
        return self._data[:, 1]

    @classmethod
    def from_output(cls, filename):
        """Extract native two-column records from a SPLEED *_SPEC.out log.

        Skip unrelated diagnostic lines; do not interpret seven-column SPC
        records here. The SPLOUT diagram marker precedes the first result.
        Return None when no two-column records are present.
        """
        rows = []
        in_diagram = False
        with open(filename, encoding="utf-8", errors="replace") as output:
            for line in output:
                if _DIAGRAM.search(line):
                    in_diagram = True
                    continue
                if not in_diagram:
                    continue
                parts = line.split()
                if len(parts) != 2 or not all(_D_NUMBER.fullmatch(p) for p in parts):
                    continue
                rows.append([float(p.replace("d", "e").replace("D", "E")) for p in parts])
        return cls(rows) if rows else None
