from io import StringIO

import numpy as np

from ..definitions.sfn import SFNOutputFile, definition
from ..output_files import OutputFile


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
    1    2
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


def test_sfn_grammar_and_generated_mesh_view():
    output = definition.read_from_string(SFN, alat=2.0)

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
    assert first.rmt == 0.6
    assert first.faces[0].NVERT() == 3
    assert first.faces[0].VERTICES().shape == (3, 3)

    second = output.mesh_for_idx(2)
    assert second is not None
    assert second.sfn.shape == (1, 2)
    assert output.mesh_for_idx(3) is None


def test_sfn_write_roundtrip():
    output = definition.read_from_string(SFN)
    written = output.to_string(validate=True)
    reparsed = definition.read_from_file(StringIO(written))

    assert reparsed.NM() == output.NM()
    for original, copy in zip(output.meshes, reparsed.meshes):
        np.testing.assert_array_equal(copy.jrcut, original.jrcut)
        np.testing.assert_allclose(copy.radial_mesh, original.radial_mesh)
        np.testing.assert_array_equal(copy.sfn_lm, original.sfn_lm)
        np.testing.assert_allclose(copy.sfn, original.sfn)
        np.testing.assert_allclose(
            copy.faces[0].VERTICES(), original.faces[0].VERTICES()
        )


def test_sfn_is_registered(tmp_path):
    path = tmp_path / "shape.sfn"
    path.write_text(SFN)

    output = OutputFile.from_file(path, unknown=False)

    assert isinstance(output, SFNOutputFile)
    assert output.NM() == 2
