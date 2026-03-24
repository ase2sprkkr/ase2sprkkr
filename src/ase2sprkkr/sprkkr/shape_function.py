"""
Reader for SPR-KKR shape function files (``*.sfn``).

The shape function file is a pure Fortran columnar file — no keyword=value
headers — so a custom line-oriented reader is used rather than the ase2sprkkr
grammar engine.

Public API
----------
``read_shape_function(path, alat=1.0)``
    Read an sfn file and return a :class:`ShapeFunction` object.

``ShapeFunction``
    Container for all NM mesh blocks.  Indexable: ``sfn[im]`` returns the
    :class:`ShapeFunctionMesh` for mesh *im* (0-based).

``ShapeFunctionMesh``
    Data for one mesh block:

    - ``idx``        – 1-based mesh index
    - ``npan``       – number of integration panels
    - ``nr``         – number of radial points
    - ``jrcut``      – panel boundary indices (length NPAN)
    - ``rmesh``      – radial coordinates in Bohr (length NR)
    - ``nsfn``       – number of shape-function components
    - ``sfn_lm``     – LM index for each component (length NSFN)
    - ``sfn``        – shape function values, shape (NSFN, NR)
    - ``rmt``        – muffin-tin radius [Bohr]
    - ``rmtfill``    – muffin-tin fill radius [alat]
    - ``vol``        – Voronoi cell volume [alat^3]
    - ``nface``      – number of Voronoi faces
    - ``faces``      – list of face dicts: ``{'nvert': int, 'verts': ndarray(NVERT,3), 'acoef': ndarray(4)}``
"""

import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class ShapeFunctionMesh:
    """Shape function data for one radial mesh block in the SFN file."""

    idx:     int
    npan:    int                             # number of integration panels
    nr:      int                             # number of SFN radial mesh points

    jrcut:   List[int]      = field(default_factory=list)
    rmesh:   np.ndarray     = field(default_factory=lambda: np.array([]))
    nsfn:    int            = 0
    sfn_lm:  List[int]      = field(default_factory=list)
    sfn:     Optional[np.ndarray] = None    # shape (NSFN, NR)

    # Summary section (filled after all mesh data)
    rmt:     float = 0.   # RMTRED — muffin-tin radius [Bohr]
    rmtfill: float = 0.   # RMTFILL [alat]
    vol:     float = 0.   # Voronoi cell volume [alat^3]
    nface:   int   = 0
    faces:   list  = field(default_factory=list)
    # Each face: {'nvert': int, 'verts': ndarray(NVERT,3), 'acoef': ndarray(4)}

    def __repr__(self):
        return (f"ShapeFunctionMesh(idx={self.idx}, npan={self.npan}, "
                f"nr={self.nr}, nsfn={self.nsfn}, nface={self.nface})")

    # ── spherical-harmonic reconstruction ────────────────────────────────────

    @staticmethod
    def _index2lm(index: int):
        """Convert SPR-KKR combined lm index (1-based) to (l, m)."""
        n = -1; L = 0; found = False
        while not found:
            n += 1
            if n ** 2 + 1 > index:
                L = n - 1; found = True
        return L, -L + (index - L ** 2 - 1)

    @staticmethod
    def _real_sph_harm(m: int, l: int, phi, theta):
        """Real spherical harmonic Y_l^m in the SPR-KKR convention."""
        from scipy.special import sph_harm as _sph_harm
        if m == 0:
            return np.real(_sph_harm(0, l, phi, theta))
        elif m > 0:
            return np.real(
                1.0 / np.sqrt(2) * (
                    _sph_harm(-m, l, phi, theta) + (-1) ** m * _sph_harm(m, l, phi, theta)
                )
            )
        else:
            return np.real(
                1j / np.sqrt(2) * (
                    _sph_harm(-m, l, phi, theta) - (-1) ** m * _sph_harm(m, l, phi, theta)
                )
            )

    def to_3d_grid(self, n: int = 80):
        """
        Reconstruct the shape function on an (n×n×n) Cartesian grid by
        summing all lm components with real spherical harmonics::

            Θ(r) = Σ_lm  sfn_lm(r) × Y_lm_real(θ, φ)

        The grid spans ±rmax in all directions where rmax = rmesh.max() * 1.02.

        Returns
        -------
        xyz : ndarray shape (n, 3) — grid coordinates along one axis [Bohr]
        grid : ndarray shape (n, n, n) — Θ values
        """
        rmax = self.rmesh[-1] * 1.02
        lin  = np.linspace(-rmax, rmax, n)
        X, Y, Z = np.meshgrid(lin, lin, lin, indexing='ij')

        R     = np.sqrt(X**2 + Y**2 + Z**2)
        # avoid division by zero at origin
        with np.errstate(invalid='ignore', divide='ignore'):
            theta = np.where(R > 0, np.arccos(np.clip(Z / R, -1, 1)), 0.)
            phi   = np.where(R > 0, np.sign(Y) * np.arccos(
                        np.clip(X / np.sqrt(np.maximum(X**2 + Y**2, 1e-30)), -1, 1)), 0.)

        grid = np.zeros((n, n, n), dtype=float)
        for i, lm_idx in enumerate(self.sfn_lm):
            l, m = self._index2lm(lm_idx)
            radial = np.interp(R, self.rmesh, self.sfn[i], left=0., right=0.)
            grid  += radial * self._real_sph_harm(m, l, phi, theta)

        return lin, grid


class ShapeFunction:
    """
    Complete parsed shape function (``*.sfn``) file.

    Parameters
    ----------
    meshes : list of ShapeFunctionMesh
    alat   : lattice constant [Bohr] used during parsing
    """

    def __init__(self, meshes: List[ShapeFunctionMesh], alat: float = 1.0):
        self._meshes = meshes
        self.alat    = alat

    @property
    def nm(self) -> int:
        """Number of mesh blocks."""
        return len(self._meshes)

    def __len__(self):
        return len(self._meshes)

    def __getitem__(self, idx):
        return self._meshes[idx]

    def __iter__(self):
        return iter(self._meshes)

    def mesh_for_idx(self, idx: int) -> Optional[ShapeFunctionMesh]:
        """Return the mesh with 1-based index *idx*, or None."""
        for m in self._meshes:
            if m.idx == idx:
                return m
        return None

    def __repr__(self):
        return f"ShapeFunction(nm={self.nm}, alat={self.alat:.6g} Bohr)"


# ── internal helpers ──────────────────────────────────────────────────────────

def _read_floats(lines, ptr, count):
    """Read *count* float values starting at *ptr*, consuming as many lines as
    needed (4 values per line is the SPR-KKR convention for sfn data)."""
    vals = []
    while len(vals) < count and ptr < len(lines):
        for v in lines[ptr].split():
            vals.append(float(v))
            if len(vals) == count:
                break
        ptr += 1
    return vals, ptr


# ── public reader ─────────────────────────────────────────────────────────────

def read_shape_function(path: str, alat: float = 1.0) -> ShapeFunction:
    """
    Read a SPR-KKR shape function file.

    File layout
    -----------
    Line 0       : NM  <unused>
    Per mesh k = 0 .. NM-1:
      Line       : NPAN  NR
      Lines      : panel boundary indices (ceil(NPAN/16) lines, 5 chars/index)
      Lines      : NR r-values as (r dr) pairs (ceil(NR/2) lines)
                   → r (column 0 and 2) stored; scaled by *alat* to give Bohr
      Line       : NSFN
      Per sfn component i = 0 .. NSFN-1:
        Line     : LM_index
        Lines    : SFN values (ceil(NR/4) lines, up to 4 values/line)
    Summary section (after all NM meshes, repeated for each mesh):
      Header     : "IM RMTRED0_M(IM) RMTFILL_M(IM) VOL_M(IM)"
      Line       : IM  RMTRED  RMTFILL  VOL  ...
      Header     : "NFACE_M(IM)"
      Line       : NFACE
      Per face:
        Header   : "IFC NVERT_FCM(IFC,IM)"
        Line     : IFC  NVERT
        Lines    : NVERT × 3 float vertex coordinates
        Header   : "ACOEF_FM"
        Line     : 4 floats (plane coefficients)

    Parameters
    ----------
    path : str
        Path to the ``*.sfn`` file.
    alat : float
        Lattice constant in Bohr; r-mesh values in the file are in units of
        *alat* and are converted to Bohr on reading.

    Returns
    -------
    ShapeFunction
    """
    with open(path) as fh:
        lines = fh.readlines()

    if not lines:
        return ShapeFunction([], alat)

    ptr = 0

    # ── header: NM ────────────────────────────────────────────────────────────
    nm   = int(lines[ptr].split()[0])
    ptr += 1

    meshes: List[ShapeFunctionMesh] = []

    # ── mesh data blocks ──────────────────────────────────────────────────────
    for k in range(nm):
        if ptr >= len(lines):
            break

        # NPAN  NR
        parts = lines[ptr].split(); ptr += 1
        if len(parts) < 2:
            continue
        npan = int(parts[0])
        nr   = int(parts[1])
        mesh = ShapeFunctionMesh(idx=k + 1, npan=npan, nr=nr)

        # Panel boundary indices: ceil(NPAN/16) lines, 5 chars per index
        nlines_pan = math.ceil(npan / 16)
        for _ in range(nlines_pan):
            if ptr >= len(lines):
                break
            raw = lines[ptr].rstrip('\n'); ptr += 1
            # Fixed 5-char fields
            nvals = max(0, len(raw) // 5)
            mesh.jrcut.extend(int(raw[i*5:(i+1)*5]) for i in range(nvals) if raw[i*5:(i+1)*5].strip())

        # Radial mesh: ceil(NR/2) lines, each line has 2 (r  dr) pairs
        # Take columns 0 and 2 (skip the dr values), scale by alat → Bohr
        sfn_r = np.zeros(nr)
        j = 0
        for _ in range(math.ceil(nr / 2)):
            if ptr >= len(lines) or j >= nr:
                break
            cols = lines[ptr].split(); ptr += 1
            if len(cols) >= 3:
                sfn_r[j] = float(cols[0]); j += 1
                if j < nr:
                    sfn_r[j] = float(cols[2]); j += 1
            elif len(cols) >= 1:
                sfn_r[j] = float(cols[0]); j += 1
        mesh.rmesh = sfn_r * alat

        # NSFN
        if ptr >= len(lines):
            meshes.append(mesh)
            continue
        mesh.nsfn = int(lines[ptr].split()[0]); ptr += 1

        sfn_arr = np.zeros((mesh.nsfn, nr))
        for i in range(mesh.nsfn):
            if ptr >= len(lines):
                break
            mesh.sfn_lm.append(int(lines[ptr].split()[0])); ptr += 1
            col = 0
            for _ in range(math.ceil(nr / 4)):
                if ptr >= len(lines) or col >= nr:
                    break
                for v in lines[ptr].split():
                    if col >= nr:
                        break
                    sfn_arr[i, col] = float(v)
                    col += 1
                ptr += 1
        mesh.sfn = sfn_arr
        meshes.append(mesh)

    # ── summary section ───────────────────────────────────────────────────────
    # Format: per mesh, a small block starting with a header line containing
    # the column names, followed by the data line, then face data.
    # We scan line-by-line looking for known header keywords.
    mesh_iter = iter(meshes)
    cur_mesh  = next(mesh_iter, None)

    while ptr < len(lines) and cur_mesh is not None:
        line = lines[ptr].strip()

        if line.startswith('IM RMTRED'):
            # next line: IM  RMTRED  RMTFILL  VOL ...
            ptr += 1
            if ptr < len(lines):
                parts = lines[ptr].split(); ptr += 1
                if len(parts) >= 4:
                    cur_mesh.rmt     = float(parts[1]) * alat
                    cur_mesh.rmtfill = float(parts[2])
                    cur_mesh.vol     = float(parts[3])

        elif line.startswith('NFACE_M'):
            ptr += 1
            if ptr < len(lines):
                cur_mesh.nface = int(lines[ptr].split()[0]); ptr += 1

        elif line.startswith('IFC NVERT'):
            ptr += 1
            if ptr < len(lines):
                parts = lines[ptr].split(); ptr += 1
                nvert = int(parts[1]) if len(parts) >= 2 else 0
                verts = np.zeros((nvert, 3))
                for vi in range(nvert):
                    if ptr < len(lines):
                        verts[vi] = [float(x) for x in lines[ptr].split()[:3]]
                        ptr += 1
                cur_mesh.faces.append({'nvert': nvert, 'verts': verts, 'acoef': None})

        elif line.startswith('ACOEF_FM'):
            ptr += 1
            if ptr < len(lines) and cur_mesh.faces:
                cur_mesh.faces[-1]['acoef'] = np.array([float(x) for x in lines[ptr].split()[:4]])
                ptr += 1
            # After the last face of this mesh, advance to next mesh
            if cur_mesh.faces and len(cur_mesh.faces) == cur_mesh.nface:
                cur_mesh = next(mesh_iter, None)

        else:
            ptr += 1

    return ShapeFunction(meshes, alat)


__all__ = ['ShapeFunction', 'ShapeFunctionMesh', 'read_shape_function']
