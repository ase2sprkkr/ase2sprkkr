from io import StringIO
import os
import sys
from types import ModuleType

import numpy as np
import pytest
from ase.build import bulk

from ..definitions.sfn import SFNOutputFile, definition
from ..output_files import OutputFile
from ...potentials.potentials import Potential
from ...sprkkr.calculator import SPRKKR
from ...common.warnings import DataValidityError


SFN = """    2    1    3
    1    3
    3
  1.00000000000000E-01  1.00000000000000E-02  2.00000000000000E-01  1.00000000000000E-02
  3.00000000000000E-01  1.00000000000000E-02
    2
    1
  1.00000000000000E+00  2.00000000000000E+00  3.00000000000000E+00
    5
  4.00000000000000E+00  5.00000000000000E+00  6.00000000000000E+00
    2    2
    1    1
  4.00000000000000E-01  2.00000000000000E-02  5.00000000000000E-01  2.00000000000000E-02
    1
    9
  7.00000000000000E+00  8.00000000000000E+00
IM RMTRED0_M(IM) RMTFILL_M(IM) VOL_M(IM)
         1  3.00000000000000000E-01  2.00000000000000000E+00  1.00000000000000000E+01
NFACE_M(IM)
         1
IFC NVERT_FCM(IFC,IM)
         1         3
  0.00000000000000000E+00  0.00000000000000000E+00  0.00000000000000000E+00
  1.00000000000000000E+00  0.00000000000000000E+00  0.00000000000000000E+00
  0.00000000000000000E+00  1.00000000000000000E+00  0.00000000000000000E+00
ACOEF_FM
  1.00000000000000000E+00  0.00000000000000000E+00  0.00000000000000000E+00  0.00000000000000000E+00
IM RMTRED0_M(IM) RMTFILL_M(IM) VOL_M(IM)
         2  5.00000000000000000E-01  3.00000000000000000E+00  2.00000000000000000E+01
NFACE_M(IM)
         1
IFC NVERT_FCM(IFC,IM)
         1         2
  0.00000000000000000E+00  0.00000000000000000E+00  1.00000000000000000E+00
  0.00000000000000000E+00  1.00000000000000000E+00  1.00000000000000000E+00
ACOEF_FM
  0.00000000000000000E+00  0.00000000000000000E+00  1.00000000000000000E+00 -1.00000000000000000E+00"""


def synthetic_potential():
    return Potential.from_atoms(bulk("Fe", "bcc", a=2.0))


def test_sfn_grammar_and_generated_mesh_view():
    output = definition.read_from_string(SFN, potential=synthetic_potential())

    assert output.nm == 2
    assert len(output.RADIAL_MESHES) == 2
    assert len(output.VORONOI_MESHES) == 2

    first = output.meshes[0]
    assert first.idx == 1
    assert first.npan == 1
    assert first.nr == 3
    assert first.nsfn == 2
    assert first.nface == 1
    np.testing.assert_array_equal(first.jrcut, [3])
    np.testing.assert_allclose(first.rmesh, [0.2, 0.4, 0.6])
    np.testing.assert_allclose(first.drmesh, [0.02, 0.02, 0.02])
    np.testing.assert_array_equal(first.sfn_lm, [1, 5])
    np.testing.assert_allclose(first.sfn, [[1, 2, 3], [4, 5, 6]])
    assert first.rmt == pytest.approx(0.6)
    assert first.faces[0].NVERT() == 3
    assert first.faces[0].VERTICES().shape == (3, 3)

    second = output.mesh_for_idx(2)
    assert second is not None
    assert second.sfn.shape == (1, 2)
    assert output.mesh_for_idx(3) is None


def test_sfn_write_roundtrip(tmp_path):
    output = definition.read_from_string(SFN)
    filename = tmp_path / "roundtrip.sfn"
    output.save_to_file(filename, validate=True)
    reparsed = OutputFile.from_file(filename, unknown=False)

    assert reparsed.NM() == output.NM()
    for original, copy in zip(output.meshes, reparsed.meshes):
        np.testing.assert_array_equal(copy.jrcut, original.jrcut)
        np.testing.assert_allclose(copy.radial_mesh, original.radial_mesh)
        np.testing.assert_array_equal(copy.sfn_lm, original.sfn_lm)
        np.testing.assert_allclose(copy.sfn, original.sfn)
        np.testing.assert_allclose(
            copy.faces[0].VERTICES(), original.faces[0].VERTICES()
        )


@pytest.mark.slow
@pytest.mark.skipif(
    os.environ.get("DO_NOT_RUN_SPRKKR", "") != "",
    reason="The test requires a working SPR-KKR executable",
)
def test_written_sfn_is_accepted_by_sprkkr(tmp_path):
    """SPR-KKR can reuse an SFN file serialized by this implementation."""

    options = {
        "NL": 3,
        "NE": 5,
        "NKTAB": 5,
        "NITER": 1,
        "FULLPOT": True,
        "MODE": "SREL",
    }
    first = SPRKKR(atoms=bulk("Al", "fcc", a=4.0)).calculate(
        directory=tmp_path,
        mpi=False,
        options=options,
        empty_spheres=False,
        print_output=False,
    )
    assert first.sfn_generated is True

    second_directory = tmp_path / "reuse"
    second_directory.mkdir()
    rewritten = second_directory / "rewritten.sfn"
    first.sfn.save_to_file(rewritten, validate=True)

    second = SPRKKR(potential=first.potential_filename).calculate(
        directory=second_directory,
        mpi=False,
        options={**options, "SFNFIL": rewritten.name},
        empty_spheres=False,
        print_output=False,
    )

    assert second.sfn_generated is False
    assert second.sfn_filename == str(rewritten)
    assert len(second.iterations) == 1


def test_sfn_dependent_array_lengths_are_validated():
    output = definition.read_from_string(SFN)

    invalid_values = (
        (output.RADIAL_MESHES[0].RMESH, np.ones((2, 2))),
        (output.RADIAL_MESHES[0].SHAPE_FUNCTIONS[0].VALUES, np.ones(2)),
        (output.VORONOI_MESHES[0].FACES[0].VERTICES, np.ones((2, 3))),
    )
    for option, value in invalid_values:
        with pytest.raises(DataValidityError, match=r"expected \(3[,)]"):
            option.set(value)


def test_sfn_is_registered(tmp_path):
    path = tmp_path / "shape.sfn"
    path.write_text(SFN)

    output = OutputFile.from_file(path, unknown=False)

    assert isinstance(output, SFNOutputFile)
    assert output.NM() == 2


def test_output_file_read_from_file_accepts_potential():
    potential = synthetic_potential()
    output = SFNOutputFile(definition=definition)

    output.read_from_file(StringIO(SFN), potential=potential)

    assert output.potential is potential
    assert output.alat == pytest.approx(2.0)


def test_sfn_shape_reconstruction_and_plot_dispatch(tmp_path):
    import inspect
    import matplotlib.pyplot as plt

    output = definition.read_from_string(SFN, potential=synthetic_potential())
    mesh = output.meshes[0]

    linear, grid = mesh.to_2d_grid(n=21)
    assert linear.shape == (21,)
    assert grid.shape == (21, 21)
    np.testing.assert_allclose(grid[10, 10], 1 / np.sqrt(4 * np.pi))
    np.testing.assert_allclose(mesh.panel_boundaries, [0.6])

    figure = plt.figure()
    mesh_axis = figure.add_subplot(131, projection="3d")
    shape_axis = figure.add_subplot(132)
    radial_axis = figure.add_subplot(133)
    assert mesh.plot(what="mesh", axis=mesh_axis) is mesh_axis
    assert mesh.plot(
        what="slice", axis=shape_axis, plane="xy", n=21, colorbar=False
    ) is shape_axis
    assert mesh.plot(what="radial", axis=radial_axis) is radial_axis
    assert mesh_axis.collections
    assert shape_axis.collections
    assert len(radial_axis.lines) >= mesh.nsfn
    plt.close(figure)

    figure = output.plot(
        what="slice",
        mesh=0,
        plane="xy",
        n=21,
        colorbar=False,
        show=False,
    )
    assert len(figure.axes) == 1
    plt.close(figure)

    common_parameters = {
        "mesh",
        "layout",
        "separate_plots",
        "figsize",
        "filename",
        "show",
        "dpi",
        "latex",
    }
    for method in (
        output.plot_mesh,
        output.plot_cell,
        output.plot_periodic,
        output.plot_shape,
        output.plot_slice,
        output.plot_radial,
    ):
        parameters = inspect.signature(method).parameters
        assert common_parameters <= parameters.keys()
        assert "what" not in parameters

    output.plot_slice(
        separate_plots=True,
        filename=tmp_path / "shape.png",
        plane="xy",
        n=21,
        colorbar=False,
        show=False,
    )
    assert (tmp_path / "shape_IM1.png").is_file()
    assert (tmp_path / "shape_IM2.png").is_file()


def test_sfn_potential_cell_and_periodic_plots(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    path = tmp_path / "Fe_SCF.sfn"
    path.write_text(SFN)
    potential = Potential.from_atoms(bulk("Fe"))
    output = OutputFile.from_file(
        path, try_only="sfn", potential=potential
    )

    assert output.potential is potential
    assert output.alat == pytest.approx(2.87)
    figure = output.plot_cell(mesh=0, show=False)
    axis = figure.axes[0]
    assert len(axis.lines) == 12
    assert axis.collections
    plt.close(figure)

    colors = output._periodic_mesh_colors(
        {1: None, 2: None, 3: None, 4: None}, None
    )
    assert len({tuple(color) for color in colors.values()}) == 4
    assert output._periodic_mesh_colors({1: None}, "red") == {1: "red"}

    figure = output.plot_periodic(mesh=0, show=False)
    axis = figure.axes[0]
    assert len(axis.lines) == 12
    assert axis.collections
    assert [text.get_text() for text in axis.get_legend().texts] == ["IMQ 1: Fe"]
    plt.close(figure)

    atoms, _, meshes, cell, inverse_cell = output._cell_context(
        (output.meshes[0],), True
    )
    assert len(atoms) > 1
    mesh_indices = np.resize(np.asarray((1, 2)), len(atoms))
    monkeypatch.setattr(
        output,
        "_cell_context",
        lambda selected, conventional: (
            atoms,
            mesh_indices,
            meshes,
            cell,
            inverse_cell,
        ),
    )
    figure = output.plot_periodic(mesh=0, show=False)
    plt.close(figure)

    without_potential = OutputFile.from_file(path, try_only="sfn")
    assert repr(without_potential) == "SFNOutputFile(nm=2, alat=unknown)"
    with pytest.raises(ValueError, match="potential is required"):
        without_potential.plot_cell(mesh=0, show=False)


def test_sfn_isosurface_in_the_cell(monkeypatch):
    import matplotlib.pyplot as plt

    measure = ModuleType("skimage.measure")
    requested_levels = []

    def marching_cubes(values, level, spacing):
        requested_levels.append(level)
        vertices = np.asarray(
            ((1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2)),
            dtype=float,
        )
        vertices *= spacing
        faces = np.asarray(
            ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3)),
            dtype=int,
        )
        return vertices, faces, None, None

    measure.marching_cubes = marching_cubes
    skimage = ModuleType("skimage")
    skimage.__path__ = []
    skimage.measure = measure
    monkeypatch.setitem(sys.modules, "skimage", skimage)
    monkeypatch.setitem(sys.modules, "skimage.measure", measure)

    potential = Potential.from_atoms(bulk("Fe"))
    output = definition.read_from_string(SFN, potential=potential)
    figure = output.plot_shape(mesh=0, n=9, show=False)
    axis = figure.axes[0]
    assert len(axis.lines) == 12
    assert requested_levels == [pytest.approx(0.5)]
    assert len(axis.collections) >= 3  # isosurface, mesh edges, and atoms
    plt.close(figure)


def test_sfn_mesh_overlay_accepts_faces_with_different_vertex_counts():
    import matplotlib.pyplot as plt

    figure = plt.figure()
    axis = figure.add_subplot(111, projection="3d")
    triangle = np.asarray(((0, 0, 0), (1, 0, 0), (0, 1, 0)))
    quadrilateral = np.asarray(
        ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    )

    SFNOutputFile._draw_polygon_edges(
        axis, (triangle, quadrilateral), "black", 0.8
    )

    assert len(axis.collections) == 2
    plt.close(figure)
