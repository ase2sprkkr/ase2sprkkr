import numpy as np

from ..radial_meshes import ExponentialMesh


def test_mesh_interpolate_maps_source_coordinates_to_mesh_coordinates():
    mesh = ExponentialMesh(1.0, np.log(2.0), 0, 0.0, 3, 0.0)
    source_coors = np.array([1.0, 2.0, 4.0, 8.0])
    values = 3.0 * source_coors + 2.0

    np.testing.assert_allclose(
        mesh.interpolate(source_coors, values),
        3.0 * mesh.coors + 2.0,
    )
