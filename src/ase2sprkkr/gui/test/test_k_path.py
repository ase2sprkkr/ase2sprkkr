"""Tests for Brillouin-zone construction used by the interactive k-path tool."""

import numpy as np
import pytest
import matplotlib.pyplot as plt
from ase import Atoms
from ase.cell import Cell

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..k_path import brillouin_zone, k_path_gui


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


@pytest.mark.parametrize(
    "cell",
    [Cell(np.eye(3) * 2.8), Cell.fromcellpar([2.8, 3.1, 3.7, 70, 80, 75])],
)
def test_k_path_gui_normalizes_cell_without_modifying_atoms(cell, monkeypatch):
    atoms = Atoms("Fe", cell=cell, pbc=True)
    original = atoms.cell.array.copy()
    figures = set(plt.get_fignums())
    drawn = []

    def show():
        # Exercise the actual plot construction and rendering without waiting
        # for interactive input. Closing without a selection returns None.
        plt.gcf().canvas.draw()
        drawn.append(True)

    monkeypatch.setattr(plt, "show", show)
    try:
        assert k_path_gui(atoms) is None
        assert drawn == [True]
        np.testing.assert_array_equal(atoms.cell.array, original)
    finally:
        for number in set(plt.get_fignums()) - figures:
            plt.close(number)
