import os
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
from ..task_result import TaskResult
from ...common.file_utils import filename_from_file
from ...common.options import Option
from ..readers.scf import ScfOutputParser, ScfResult, atomic_types_definition
from ..readers.jxc import JxcOutputReader

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)


class TestOutput(TestCase):
    def test_task_result_normalizes_existing_files(self):
        with TemporaryDirectory() as directory:
            filename = os.path.join(directory, "result.out")
            with open(filename, "wb") as output:
                output.write(b"output")
            with open(filename, "rb") as output:
                assert filename_from_file(output) == filename
                assert filename_from_file(Path(filename)) == filename
                result = TaskResult(None, None, directory, output_file=output)
                assert result.path_to("output") == filename
                assert isinstance(result.files["output"], Option)
                assert result.files["output"].path() == filename
                assert "open" in result.files["output"].actions()

            relative = TaskResult(None, None, directory, output_file="result.out")
            assert relative.path_to("output") == filename

            inferred_directory = TaskResult(None, None, None, output_file=filename)
            assert inferred_directory.directory == directory
            assert inferred_directory.path_to("output") == filename

            nested = TaskResult(None, None, None, output_file=os.path.join("subdir", "result.out"))
            assert nested.path_to("output") == os.path.join(os.getcwd(), "subdir", "result.out")

            stream = BytesIO(b"output")
            with self.assertRaises(TypeError):
                filename_from_file(stream)
            assert filename_from_file(stream, None) is None
            stream_result = TaskResult(None, None, directory, output_file=stream)
            assert stream_result.output_file is None
            assert "output" not in stream_result.files

            no_output = TaskResult(None, None, directory)
            assert no_output.output_file is None
            assert "output" not in no_output.files

    def test_scf(self):
        atomic_types_definition.parse(
            """  33 E= 0.6083 0.0000          IT=   1  Li_1
         DOS      NOS     P_spin   m_spin    P_orb    m_orb    B_val      B_core
  s    0.4387   0.0296    0.0000   0.0000   0.00000  0.00000    0.00 s      0.00
  p    1.2579   0.0962    0.0000   0.0000   0.00000  0.00000    0.00 ns     0.00
  d    0.6886   0.0926    0.0000   0.0000   0.00000  0.00000    0.00 cor    0.00
  f    0.3427   0.0476    0.0000   0.0000   0.00000  0.00000    0.00
 sum   2.7279   0.2660    0.0000   0.0000   0.00000  0.00000    0.00 v+c    0.00
 E_band         0.11559127 [Ry]
dipole moment   1      0.0000000000000000      0.0000000000000000      0.0000000000000000"""
        )  # NOQA: E122

    def test_output(self):
        path = os.path.join(os.path.dirname(__file__), "..", "examples", "scf.out")

        # read_from_file is both method and class_method
        for reader in ScfOutputParser, ScfOutputParser():
            out = ScfResult(None, None, None)
            reader.read_from_file(path, read_args=[out])

            self.assertEqual(out.iterations[-1]["energy"]["EMIN"](), -0.5)
            self.assertEqual(len(out.iterations[-1]["atomic_types"]), 2)
            self.assertEqual(out.iterations[0].converged(), False)
            self.assertEqual(out.iterations[-1].converged(), True)
            values = out.output_values
            assert isinstance(values["Converged"], Option)
            assert values["Converged"]._definition.is_generated
            self.assertEqual(values["Number of iterations"](), len(out.iterations))
            self.assertEqual(values["Fermi energy"].actions(), ("plot",))
            self.assertEqual(values["Data"].actions(), ("data",))
            assert values["Data"].data() is out

    def test_jxc(self):
        path = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", "examples"))
        result = JxcOutputReader.read_from_file("Fe_jxc.out", directory=path)
        assert result.jxc_filename == os.path.join(path, "Fe_JXC_XCPLTEN_Jij.dat")
        assert result.dij_filename == os.path.join(path, "Fe_JXC_XCPLTEN_Dij.dat")
        assert result.dmi_filename == os.path.join(path, "Fe_JXC_DMIVEC_Dij.dat")
        assert isinstance(result.output_values["jxc"], Option)
        assert result.output_values["jxc"]() == "Fe_JXC_XCPLTEN_Jij.dat"
        assert callable(result.output_values["jxc"].parsed_file)
        assert result.mean_field_curie_temperature == 763.4

    def test_reader_keeps_underlying_output_filename(self):
        directory = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", "examples"))
        filename = os.path.join(directory, "Fe_jxc.out")
        with open(filename, "rb") as output:
            result = JxcOutputReader.read_from_file(output, directory=directory)
        assert result.path_to("output") == filename
        assert result.files["output"].path() == filename
