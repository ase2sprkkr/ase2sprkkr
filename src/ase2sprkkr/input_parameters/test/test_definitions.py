import os
import re
import tempfile
from types import SimpleNamespace

import numpy as np
import pytest

from ...common.grammar_types import Keyword
from ...common.warnings import DataValidityError, DataValidityWarning
from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition


class TestDefinitions:
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
        parameters = SimpleNamespace(CONTROL=SimpleNamespace(DATASET=lambda: "Fe"))

        class DynamicDefinition(InputParametersDefinition):
            def output_files(self, parameters, directory):
                return (
                    *super().output_files(parameters, directory),
                    os.path.join(directory, "dynamic.out"),
                )

        definition = DynamicDefinition(
            "test",
            [],
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
        definition = parameters._definition
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

    def test_bsf_definitions_are_shared(self):
        definition = InputParameters.definition("BSF")
        assert definition is InputParameters.definition("BSFEK")
        assert definition is InputParameters.definition("BSFKK")

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
            np.testing.assert_equal(ip.to_dict(), parsed.to_dict())

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
        np.testing.assert_equal(bsfkk.to_dict(), parsed.to_dict())

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
        assert {value for value, _ in path_type.items()} == set(paths)
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

    def test_defaults(self):
        for name in InputParameters.definitions:
            ip = InputParameters.create_input_parameters(name)
            ip.CONTROL.POTFIL = "xxx"
            definition = ip._definition
            try:
                out = ip.to_string(validate=True)
            except DataValidityError:
                if name == "BSFEK":
                    ip.TASK.KPATH = 1
                    out = ip.to_string(validate=True)
                    with pytest.raises(DataValidityError):
                        ip.TASK.NKDIR = 2
                    ip.TASK.KPATH = None
                    ip.TASK.NKDIR = 2
                    out = ip.to_string(validate=True)
                elif name in ("BSF", "BSFKK"):
                    self.set_bsfkk_grid(ip)
                    out = ip.to_string(validate=True)
                else:
                    raise
            else:
                if name == "BSFEK":
                    pytest.fail("BSFEK should not be runnable using only default arguments")

            ip2 = definition.read_from_string(out)
            np.testing.assert_equal(ip.to_dict(), ip2.to_dict())

            if name == "SCF":
                ip.MODE.MDIR[1] = 1.0, 1.0, 1.0
                ip.MODE.MDIR[4] = 1.0, 1.0, 1.0
                ip2 = definition.read_from_string(ip.to_string())
                np.testing.assert_equal(ip.to_dict(), ip2.to_dict())

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
    def test_definitions(self, monkeypatch):
        path = os.path.join(os.path.dirname(__file__), "../examples")

        def check(parameters, name):
            assert parameters.name == name
            uname = name.upper()
            assert parameters.__class__ is InputParameters
            assert parameters["CONTROL"]["ADSI"]() == uname
            assert parameters.TASK.TASK() == uname

        for filename_entry in os.listdir(path):
            if not filename_entry.endswith(".in"):
                continue
            filename = os.path.join(path, filename_entry)
            assert filename_entry[:-3].upper() in InputParameters.definitions
            definition = InputParameters.definitions[filename_entry[:-3].upper()]
            parameters = definition.read_from_file(filename)

            name = filename_entry[:-3]
            check(parameters, name)
            parameters = InputParameters.from_file(filename)
            check(parameters, name)

            if name == "scf":
                assert "MODE" not in parameters.to_string()
            if name == "arpes":
                for task in ["AIPES", "SPLEED", "BAND"]:
                    parameters.TASK.TASK = task

        definition = InputParameters.definitions["PHAGEN"]
        task_definition = definition["TASK"]["TASK"]
        monkeypatch.setattr(task_definition, "default_value", "scf")
        parameters = definition.read_from_file(os.path.join(path, "phagen.in"))
        assert parameters.TASK.TASK() == "PHAGEN"

        with tempfile.TemporaryFile("w+") as fp:
            parameters.save_to_file(fp)
            fp.seek(0)
            out = fp.read()
            assert "TASKPHAGEN" in re.sub(r"\s", "", out)
