"""Tests for Brillouin-zone construction used by the interactive k-path tool."""

import numpy as np
import pytest
from ase.cell import Cell

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..k_path import brillouin_zone


@pytest.mark.parametrize(
    "cell, expected_faces, expected_vertices",
    [
        (Cell(np.eye(3)), 6, 8),
        (Cell.fromcellpar([1, 1, 1, 60, 60, 60]), 14, 24),
        (Cell.fromcellpar([1, 1.3, 1.7, 70, 80, 75]), 14, 24),
    ],
)
def test_brillouin_zone_topology(cell, expected_faces, expected_vertices):
    reciprocal, faces = brillouin_zone(cell)
    vertices = np.unique(np.round(np.concatenate(faces), decimals=12), axis=0)

    assert np.allclose(reciprocal, cell.reciprocal())
    assert len(faces) == expected_faces
    assert len(vertices) == expected_vertices
    # Every BZ face must be a polygon, not an open or degenerate segment.
    assert all(face.shape[1] == 3 and len(face) >= 3 for face in faces)
