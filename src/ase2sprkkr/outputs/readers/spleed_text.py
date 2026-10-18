"""Extract the SPLEED spectrum embedded in native KKRSPEC NOUT1 output.

Unlike IFILSPECOU1, NOUT1 is the calculation's _SPEC.out text log.
SPLOUT writes energy/intensity as two D12.5 fields for TYP=0 and for
non-relativistic TYP=1/2. Relativistic serial TYP=1/2 writes seven E12.5
fields. These are not standalone structured ResultFile formats.
"""

import re

import numpy as np


_DIAGRAM = re.compile(
    r"\\b(?:i\\s*\\(e\\)\\s*-?diagram|rotation diagram|scattering-angle diagram)",
    re.IGNORECASE,
)
_NUMBER = re.compile(
    r"[+-]?(?:\\d+\\.\\d*|\\.\\d+)[de][+-]?\\d{2,3}\\Z",
    re.IGNORECASE,
)


class SpleedTextSpectrum:
    """Data printed by SPLOUT inside the normal KKRSPEC output log.

    Columns: two-field rows are ENERGY, INTENSITY. Seven-field rows also
    contain THETA, ASYMMETRY, P_X, P_Y and P_Z in native SPLOUT order.
    """

    def __init__(self, data, source):
        self.data = np.asarray(data, dtype=float)
        self.source = source
        if self.data.ndim != 2 or self.data.shape[1] not in (2, 7):
            raise ValueError("SPLEED text spectrum must have two or seven columns")

    def RAW_DATA(self):
        return self.data

    def ENERGY(self):
        return self.data[:, 0]

    def INTENSITY(self):
        return self.data[:, 2 if self.data.shape[1] == 7 else 1]

    def _seven_column(self, column, name):
        if self.data.shape[1] != 7:
            raise AttributeError(f"{name} is not present in two-column SPLEED output")
        return self.data[:, column]

    def THETA(self):
        return self._seven_column(1, "THETA")

    def ASYMMETRY(self):
        return self._seven_column(3, "ASYMMETRY")

    def P_X(self):
        return self._seven_column(4, "P_X")

    def P_Y(self):
        return self._seven_column(5, "P_Y")

    def P_Z(self):
        return self._seven_column(6, "P_Z")

    @classmethod
    def from_file(cls, filename):
        """Parse the SPLOUT numerical rows; return None when absent.

        Require the SPLOUT diagram marker and Fortran exponent formats to
        avoid confusing unrelated numerical diagnostics in _SPEC.out with
        spectrum data. Other log lines are ignored.
        """
        with open(filename, encoding="utf-8", errors="replace") as stream:
            lines = stream.readlines()

        active = False
        rows = {}
        for line in lines:
            if _DIAGRAM.search(line):
                active = True
                continue
            if not active:
                continue
            tokens = line.split()
            if len(tokens) not in (2, 7):
                continue
            if not all(_NUMBER.fullmatch(value) for value in tokens):
                continue
            if len(tokens) == 2 and not all("D" in value.upper() for value in tokens):
                continue
            if len(tokens) == 7 and not all("E" in value.upper() for value in tokens):
                continue
            rows.setdefault(len(tokens), []).append(
                [float(value.replace("d", "e").replace("D", "E")) for value in tokens]
            )

        if not rows:
            return None
        if len(rows) != 1:
            raise ValueError("Mixed two- and seven-column SPLEED records in NOUT1")
        return cls(next(iter(rows.values())), filename)
