"""
Parser for *_SPEC.out files produced by the KKRSPEC photoemission code.
"""

import re
import numpy as np
from typing import Dict, Optional

from ..task_result import TaskResult, KkrOutputReader
from .default import DefaultOutputParser
from ...common.decorators import cached_property
from ...output_files.output_files import OutputFile

# Scientific-notation float, including negative numbers glued without space
# e.g. "0.0000000E+00-0.1000000E+01" (SPR-KKR omits the separator)
_FORT = r'[-+]?\d+\.\d+[Ee][+-]\d+'


class SpecResult(TaskResult):
    """
    Parsed results from a ``*_SPEC.out`` file (KKRSPEC photoemission output).

    All properties are lazily parsed from the file content on first access.

    Attributes
    ----------
    potential_barrier : dict
        Raw dict with keys ``ibar``, ``epsx``, ``zparup``, ``zpardn``,
        ``bparp`` (the last three as lists of three floats).
    stokes : dict
        Stokes vector of the photon field with keys ``s0``, ``s1_pct``,
        ``s2_pct``, ``s3_pct``.
    polarization_type : str
        Derived polarization label: ``'C+'``, ``'C-'``, ``'LH'``, or ``'LV'``.
    photon_wavevector : np.ndarray or None
        Vacuum wavevector of the photon field, shape (3,), units Bohr⁻¹.
    vector_potential : dict or None
        Jones vector with keys ``re`` and ``im``, each shape (3,).
    lattice_constants : dict
        Lattice constants, currently ``{'a': float}`` in Bohr.
    basis_vectors : dict
        ``{'real': np.ndarray (2,2), 'reciprocal': np.ndarray (2,2)}``.
    bulkrepeat : float
        Bulk interlayer repeat distance (alat units).
    spectral_data : np.ndarray
        Raw spectral table parsed after the Stokes-vector header line,
        shape (N, 6).  Empty array if the section is absent.
    """

    # ── file content ──────────────────────────────────────────────────────────

    def _get_file_content(self) -> str:
        if not hasattr(self, '_file_content'):
            if hasattr(self, 'output_file') and self.output_file:
                with open(self.output_file, 'r') as f:
                    self._file_content = f.read()
            else:
                self._file_content = ''
        return self._file_content

    # ── potential barrier ─────────────────────────────────────────────────────

    @cached_property
    def potential_barrier(self) -> Dict:
        """Potential barrier parameters from the SPEC.out header."""
        content = self._get_file_content()
        result = {}

        m = re.search(r'#\s*ibar:\s*(\d+)', content)
        if m:
            result['ibar'] = int(m.group(1))

        m = re.search(r'#\s*epsx:\s*([\d.E+-]+)', content)
        if m:
            result['epsx'] = float(m.group(1))

        # SPR-KKR writes negative numbers without whitespace separation,
        # e.g. "0.0000000E+00-0.1000000E+01 0.1000000E+01".  Use findall on
        # the full scientific-notation pattern to handle this correctly.
        for key in ('zparup', 'zpardn', 'bparp'):
            m = re.search(rf'#\s*{key}:\s*(.*)', content)
            if m:
                vals = re.findall(_FORT, m.group(1))
                if len(vals) >= 3:
                    result[key] = [float(v) for v in vals[:3]]

        return result

    # ── Stokes vector ─────────────────────────────────────────────────────────

    @cached_property
    def stokes(self) -> Dict:
        """
        Stokes vector of the photon field.

        Keys: ``s0`` (total intensity), ``s1_pct``, ``s2_pct``, ``s3_pct``
        (linear horizontal %, linear 45° %, circular %).
        """
        content = self._get_file_content()
        m = re.search(
            r'stookes-vector\s+of\s+the\s+photonfield\s+'
            r'([\d.E+-]+)\s+([\d.E+-]+)%\s+([\d.E+-]+)%\s+([-\d.E+-]+)%',
            content, re.IGNORECASE)
        if m:
            return {
                's0':     float(m.group(1)),
                's1_pct': float(m.group(2)),
                's2_pct': float(m.group(3)),
                's3_pct': float(m.group(4)),
            }
        return {'s0': 0., 's1_pct': 0., 's2_pct': 0., 's3_pct': 0.}

    @cached_property
    def polarization_type(self) -> str:
        """
        Derived polarization label from the Stokes vector.

        Returns ``'C+'``, ``'C-'``, ``'LH'``, or ``'LV'``.
        """
        s = self.stokes
        if abs(s['s3_pct']) > 50.:
            return 'C+' if s['s3_pct'] < 0 else 'C-'
        if abs(s['s1_pct']) > 50.:
            return 'LH'
        return 'LV'

    # ── photon field vectors ──────────────────────────────────────────────────

    @cached_property
    def photon_wavevector(self) -> Optional[np.ndarray]:
        """Vacuum wavevector of the photon field, shape (3,), units Bohr⁻¹."""
        content = self._get_file_content()
        m = re.search(
            r'vacuum-wavevector of the photonfield:\s+'
            r'([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)',
            content)
        if m:
            return np.array([float(x) for x in m.groups()])
        return None

    @cached_property
    def vector_potential(self) -> Optional[Dict]:
        """
        Jones vector of the photon field.

        Returns ``{'re': np.ndarray (3,), 'im': np.ndarray (3,)}`` or
        ``None`` if the line is absent.
        """
        content = self._get_file_content()
        m = re.search(
            r'aa\(i=1,3\)=\s*\(\s*([-\d.E+]+)\s+([-\d.E+]+)\)\s*'
            r'\(\s*([-\d.E+]+)\s+([-\d.E+]+)\)\s*'
            r'\(\s*([-\d.E+]+)\s+([-\d.E+]+)\)',
            content)
        if m:
            v = [float(x) for x in m.groups()]
            return {
                're': np.array([v[0], v[2], v[4]]),
                'im': np.array([v[1], v[3], v[5]]),
            }
        return None

    # ── lattice ───────────────────────────────────────────────────────────────

    @cached_property
    def lattice_constants(self) -> Dict:
        """Lattice constants parsed from the SPEC.out header."""
        content = self._get_file_content()
        result = {}
        m = re.search(r'lattice constant\s*:\s*([\d.]+)', content)
        if m:
            result['a'] = float(m.group(1))
        return result

    @cached_property
    def basis_vectors(self) -> Dict:
        """
        2D real and reciprocal basis vectors.

        Returns ``{'real': np.ndarray (2,2), 'reciprocal': np.ndarray (2,2)}``.
        """
        content = self._get_file_content()
        out = {'real': np.zeros((2, 2)), 'reciprocal': np.zeros((2, 2))}
        for m in re.finditer(
                r'(real|reciprocal) basis\s+(\d+)\s*:\s*([-\d.E+]+)\s+([-\d.E+]+)',
                content):
            out[m.group(1)][int(m.group(2)) - 1] = [float(m.group(3)),
                                                     float(m.group(4))]
        return out

    @cached_property
    def bulkrepeat(self) -> float:
        """Bulk interlayer repeat distance (alat units)."""
        content = self._get_file_content()
        m = re.search(r'bulkrepeat unit=\s*([\d.E+-]+)', content)
        return float(m.group(1)) if m else 0.

    # ── spectral data ─────────────────────────────────────────────────────────

    @cached_property
    def spectral_data(self) -> np.ndarray:
        """
        Raw spectral table from the SPEC.out file, shape (N, 6).

        Columns: θ, E_bind, I_total, I_up, I_down, I_pol (order as written
        by KKRSPEC after the Stokes-vector header line).
        Returns an empty array if the section is absent.
        """
        content = self._get_file_content()
        header_end = re.search(r'stookes-vector\s+of the photonfield', content)
        if not header_end:
            return np.array([])
        data = []
        for line in content[header_end.end():].split('\n'):
            if re.match(
                    r'^\s*[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+\s+'
                    r'[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E]+\s*$', line):
                data.append([float(x) for x in line.split()])
        return np.array(data) if data else np.array([])


class SpecOutputParser(DefaultOutputParser):
    """Reader for spec.out files."""

    async def read_output(self, stdout, result):
        await super().read_output(stdout, result)
        if hasattr(result, 'output_lines'):
            result.output_content = '\n'.join(result.output_lines)
        return result


class SpecOutputReader(KkrOutputReader):
    """Process class for spectral function calculations."""

    result_class = SpecResult
    parser_class = SpecOutputParser
