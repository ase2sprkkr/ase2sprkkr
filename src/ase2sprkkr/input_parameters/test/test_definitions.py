import os
import tempfile
import re
from types import SimpleNamespace
import pytest

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..input_parameters import InputParameters
    from ..input_parameters_definitions import InputParametersDefinition
    from ..task_definition import TaskDefinition
    from ...common.warnings import DataValidityError, DataValidityWarning
    from ...common.grammar_types import Keyword


class TestDefinitions(TestCase):
    @staticmethod
    def set_bsfkk_grid(ip):
        ip.TASK.NK1 = 10
        ip.TASK.K1 = [1.0, 0.0, 0.0]
        ip.TASK.NK2 = 10
        ip.TASK.K2 = [0.0, 1.0, 0.0]

    def test_set_option(self):
        ip = InputParameters.create("SCF")

        ip.set_option("DATASET", "changed")
        assert ip.CONTROL.DATASET() == "changed"

        with pytest.raises(KeyError, match="No option with name UNKNOWN_OPTION"):
            ip.set_option("UNKNOWN_OPTION", 1)

    def test_task_file_templates_and_overrides(self, tmp_path):
        parameters = SimpleNamespace(
            CONTROL=SimpleNamespace(DATASET=lambda: "Fe")
        )
        schema = InputParametersDefinition("test", [])

        class DynamicTaskDefinition(TaskDefinition):
            def output_files(self, parameters, directory):
                return (
                    *super().output_files(parameters, directory),
                    os.path.join(directory, "dynamic.out"),
                )

        definition = DynamicTaskDefinition(
            schema,
            name="TEST",
            input_file_templates=("{dataset}.aux",),
            output_file_templates=("{dataset}.dat",),
        )

        assert definition.input_files(parameters, str(tmp_path)) == (
            str(tmp_path / "Fe.aux"),
        )
        assert definition.output_files(parameters, str(tmp_path)) == (
            str(tmp_path / "Fe.dat"),
            str(tmp_path / "dynamic.out"),
        )

    def test_input_parameters_delegate_file_manifests(self, tmp_path, monkeypatch):
        parameters = InputParameters.create("DOS")
        definition = parameters.task_definition
        calls = []

        def input_files(actual_parameters, directory):
            calls.append(("input", actual_parameters, directory))
            return ("input.dat",)

        def output_files(actual_parameters, directory):
            calls.append(("output", actual_parameters, directory))
            return ("output.dat",)

        monkeypatch.setattr(definition, "input_files", input_files)
        monkeypatch.setattr(definition, "output_files", output_files)

        assert parameters.input_files(str(tmp_path)) == ("input.dat",)
        assert parameters.output_files(str(tmp_path)) == ("output.dat",)
        assert calls == [
            ("input", parameters, str(tmp_path)),
            ("output", parameters, str(tmp_path)),
        ]

    def test_task_binding_is_authoritative(self):
        parameters = InputParameters.create("SCF")

        assert parameters.task_name == "SCF"
        assert parameters._definition is parameters.task_definition.input_parameters_definition
        assert parameters.task_definition.executable == "kkrscf"
        assert parameters.task_definition.mpi is True

    def test_copy_preserves_task_binding(self):
        parameters = InputParameters.create("BSFEK")
        copied = parameters.copy(copy_values=True)

        assert copied.task_definition is parameters.task_definition
        assert copied.task_name == "BSFEK"
        assert copied._definition is copied.task_definition.input_parameters_definition
        assert copied.ENERGY.NE().tolist() == [200]

    def test_change_task_changes_task_and_schema_atomically(self):
        parameters = InputParameters.create("DOS")
        parameters.change_task("SCF")

        assert parameters.task_name == "SCF"
        assert parameters.task_definition.name == "SCF"
        assert parameters._definition is parameters.task_definition.input_parameters_definition

    def test_direct_schema_objects_are_detached(self):
        schema = InputParametersDefinition.definition_from_dict(
            {"TEST": []}
        )
        parameters = schema.create_object()

        assert parameters.task_definition is None
        assert parameters.task_name == schema.name
        with pytest.raises(RuntimeError, match="not bound"):
            parameters.is_mpi()

    def test_bsf_definitions_are_shared(self):
        definition = InputParameters.definition("BSF")
        assert definition is InputParameters.definition("BSFEK")
        assert definition is InputParameters.definition("BSFKK")

        bsf = InputParameters.task_registry["BSF"]
        bsfek = InputParameters.task_registry["BSFEK"]
        bsfkk = InputParameters.task_registry["BSFKK"]
        assert bsfek is not bsf
        assert bsfkk is not bsf
        assert bsfek.input_parameters_definition is definition
        assert bsfkk.input_parameters_definition is definition
        assert bsfek.name == "BSFEK"
        assert bsfkk.name == "BSFKK"

    def test_bsf_mode_is_inferred_when_parsing(self):
        definition = InputParameters.definition("BSFEK")
        for name in ("BSFEK", "BSFKK"):
            ip = InputParameters.create(name)
            ip.CONTROL.POTFIL = "x"
            if name == "BSFEK":
                ip.TASK.KPATH = 1
            else:
                self.set_bsfkk_grid(ip)
                ip.TASK.KA = [[0.0, 0.0, 0.0]]
            parsed = definition.read_from_string(ip.to_string(validate=True))
            assert ip.name == "bsf"
            assert parsed.name == "bsf"
            self.assertEqual(ip.to_dict(), parsed.to_dict())

    def test_bsfek_k_path_is_checked_when_saving(self):
        ip = InputParameters.create("BSFEK")
        ip.CONTROL.POTFIL = "x"

        assert ip.name == "bsf"
        assert ip.ENERGY.NE().tolist() == [200]

        with pytest.raises(DataValidityError, match="TASK.KPATH or TASK.KA"):
            ip.validate("save")

        ip.TASK.KPATH = 1
        ip.validate("save")

    def test_bsf_mode_is_selected_by_ne(self):
        bsfkk = InputParameters.create("BSFKK")

        assert bsfkk.name == "bsf"
        assert bsfkk.ENERGY.NE().tolist() == [1]
        with pytest.warns(DataValidityWarning, match="BSF-EK.*BSF-KK"):
            bsfkk.ENERGY.NE = 2
        assert bsfkk.ENERGY.NE().tolist() == [2]

        bsfek = InputParameters.create("BSFEK")
        assert bsfek.ENERGY.NE().tolist() == [200]
        with pytest.warns(DataValidityWarning, match="BSF-KK.*BSF-EK"):
            bsfek.ENERGY.NE = 1
        assert bsfek.ENERGY.NE().tolist() == [1]

    def test_bsf_modes_cannot_be_mixed(self):
        ip = InputParameters.create("BSFEK")
        ip.TASK.KPATH = 1

        with pytest.raises(DataValidityError, match="cannot be used in BSF-EK"):
            ip.TASK["K1"].set([1.0, 0.0, 0.0])

        ip = InputParameters.create("BSFKK")
        with pytest.raises(DataValidityError, match="cannot be used in BSF-KK"):
            ip.TASK["KPATH"].set(1)

    def test_bsf_mode_conflict_warns_when_parsing(self):
        definition = InputParameters.definition("BSFEK")

        bsfek = InputParameters.create("BSFEK")
        bsfek.CONTROL.POTFIL = "x"
        bsfek.TASK.KPATH = 1
        text = bsfek.to_string(validate=True).replace(
            "\tKPATH=1", "\tKPATH=1\n\tNK1=10"
        )
        with pytest.warns(DataValidityError, match="cannot be used in BSF-EK"):
            definition.read_from_string(text)

        bsfkk = InputParameters.create("BSFKK")
        bsfkk.CONTROL.POTFIL = "x"
        self.set_bsfkk_grid(bsfkk)
        text = bsfkk.to_string(validate=True).replace(
            "\tNK1=10", "\tNK=300\n\tNK1=10"
        )
        with pytest.warns(DataValidityError, match="cannot be used in BSF-KK"):
            definition.read_from_string(text)

    def test_bsfkk_grid_is_required(self):
        ip = InputParameters.create("BSFKK")
        ip.CONTROL.POTFIL = "x"

        with pytest.raises(DataValidityError, match="NK1.*NK2.*K1.*K2.*required"):
            ip.validate("save")

        self.set_bsfkk_grid(ip)
        ip.validate("save")

    def test_bsf_ka_numbering_depends_on_mode(self):
        bsfkk = InputParameters.create("BSFKK")
        bsfkk.CONTROL.POTFIL = "x"
        self.set_bsfkk_grid(bsfkk)
        bsfkk.TASK.KA = [[0.0, 0.0, 0.0]]
        out = bsfkk.to_string(validate=True)
        assert "\n\tKA={0.0,0.0,0.0}" in out
        assert "\n\tKA1=" not in out

        parsed = InputParameters.definition("BSF").read_from_string(out)
        self.assertEqual(bsfkk.to_dict(), parsed.to_dict())

        bsfek = InputParameters.create("BSFEK")
        bsfek.CONTROL.POTFIL = "x"
        bsfek.TASK.NKDIR = 1
        out = bsfek.to_string(validate=True)
        assert "\n\tKA1=" in out
        assert "\n\tKE1=" in out

    def test_bsf_path_limits(self):
        ip = InputParameters.create("BSFEK")
        ip.CONTROL.POTFIL = "x"
        path_type = ip.TASK.KPATH._definition.type
        assert isinstance(path_type, Keyword)
        paths = (0, 1, 2, 3, 4, 5, 6, 7, 10)
        assert {value for value, _ in path_type.items()} == {path for path in paths}
        for path in paths:
            for value in (path, str(path)):
                ip.TASK.KPATH = value
                assert ip.TASK.KPATH() == str(path)
                text = ip.to_string(validate=True)
                assert f"\n\tKPATH={path}\n" in text
                parsed = InputParameters.definition("BSF").read_from_string(text)
                assert parsed.TASK.KPATH() == str(path)
        for path in (8, 9, 11, "custom"):
            with pytest.raises(DataValidityError):
                ip.TASK.KPATH = path
            assert ip.TASK.KPATH() == "10"

        ip.TASK.KPATH = None
        with pytest.raises(DataValidityError, match="NKDIR cannot be greater than 9"):
            ip.TASK.NKDIR = 10

    def test_change_task(self):
        ip = InputParameters.create("DOS")
        ip.ENERGY.EMIN = 2
        ip.ENERGY.EMAX = 5
        ip.change_task("SCF", retain_values=True)
        assert ip.task_name == "SCF"
        assert ip.ENERGY.EMIN() == 2

    def test_jxc_warnings_survive_task_switch(self):
        ip = InputParameters.create("SCF")
        ip.MODE.MODE = "nrel"
        with pytest.warns(UserWarning):
            ip.change_task("jxc", retain_values=True)
        with pytest.warns(UserWarning):
            ip.MODE.MODE = "srel"
        ip.MODE.MODE = "frel"
        ip.change_task("scf")

    def test_defaults(self):
        for i in InputParameters.definitions:
            ip = InputParameters.create_input_parameters(i)
            ip.CONTROL.POTFIL = "xxx"
            df = ip._definition
            try:
                out = ip.to_string(validate=True)
            except DataValidityError:
                if i == "BSFEK":
                    ip.TASK.KPATH = 1
                    out = ip.to_string(validate=True)
                    with pytest.raises(Exception):
                        ip.TASK.NKDIR = 2
                    ip.TASK.KPATH = None
                    ip.TASK.NKDIR = 2
                    out = ip.to_string(validate=True)
                elif i in ("BSF", "BSFKK"):
                    self.set_bsfkk_grid(ip)
                    out = ip.to_string(validate=True)
                else:
                    raise
            else:
                if i == "BSFEK":
                    raise Exception("This tasks should not be runnable using the defaults argument")

            ip2 = df.read_from_string(out)
            self.assertEqual(ip.to_dict(), ip2.to_dict())

            if i == "SCF":
                ip.MODE.MDIR[1] = 1.0, 1.0, 1.0
                ip.MODE.MDIR[4] = 1.0, 1.0, 1.0
                ip2 = df.read_from_string(ip.to_string())
                self.assertEqual(ip.to_dict(), ip2.to_dict())

    def test_sfn_filename_default(self):
        ip = InputParameters.create("SCF")
        ip.CONTROL.DATASET = "calculation"
        ip.CONTROL.POTFIL = "potential.pot"

        assert ip.CONTROL.SFNFIL() == "calculation.sfn"
        assert "SFNFIL" not in ip.CONTROL.to_string()

        ip.CONTROL.DATASET = "renamed"
        assert ip.CONTROL.SFNFIL() == "renamed.sfn"
        assert "SFNFIL" not in ip.CONTROL.to_string()

        ip.CONTROL.SFNFIL = "custom.sfn"
        assert "SFNFIL=custom.sfn" in ip.CONTROL.to_string()

    @pytest.mark.slow
    def test_definitions(self):
        path = os.path.join(os.path.dirname(__file__), "../examples")

        def check(t, name):
            self.assertEqual(t.name, name)
            uname = name.upper()
            self.assertEqual(t.__class__, InputParameters)
            self.assertEqual(t["CONTROL"]["ADSI"](), uname)
            self.assertEqual(t.TASK.TASK(), uname)

        for i in os.listdir(path):
            try:
                if not i.endswith(".in"):
                    continue
                filename = os.path.join(path, i)
                self.assertTrue(i[:-3].upper() in InputParameters.definitions)
                td = InputParameters.definitions[i[:-3].upper()]
                t = td.read_from_file(filename)

                name = i[:-3]
                check(t, name)
                t = InputParameters.from_file(filename)
                check(t, name)
                assert t.task_definition is not None
                assert t._definition is t.task_definition.input_parameters_definition

                if name == "scf":
                    self.assertFalse("MODE" in t.to_string())
                if name == "arpes":
                    for task_name in ["AIPES", "SPLEED", "BAND"]:
                        t.TASK.TASK = task_name

            except Exception as e:
                raise Exception(f'Parsing of "{i}" failed with the reason: \n {e}').with_traceback(e.__traceback__)

        td = InputParameters.definitions["PHAGEN"]
        td["TASK"]["TASK"].default_value = "scf"
        t = td.read_from_file(os.path.join(path, "phagen.in"))
        self.assertEqual(t.TASK.TASK(), "PHAGEN")

        with tempfile.TemporaryFile("w+") as fp:
            t.save_to_file(fp)
            fp.seek(0)
            out = fp.read()
            assert "TASKPHAGEN" in re.sub(r"\s", "", out)
        td["TASK"]["TASK"].default_value = "PHAGEN"
