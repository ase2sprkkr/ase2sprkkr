import os
import pickle
import re
import tempfile

import numpy as np
import pytest

from ...common.grammar_types import Keyword
from ...common.warnings import DataValidityError
from ..input_parameters import InputParameters
from ..input_parameters_definitions import InputParametersDefinition
from ..task_definition import TaskDefinition


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

    def test_adsi_is_customizable_string_with_task_default(self):
        ip = InputParameters.create("SCF")

        assert ip.CONTROL.ADSI() == "SCF"

        ip.CONTROL.ADSI = "CUSTOM"
        assert ip.CONTROL.ADSI() == "CUSTOM"

        parsed = InputParameters.definition("SCF").read_from_string(ip.to_string())
        assert parsed.CONTROL.ADSI() == "CUSTOM"

    def test_task_files_and_schema_overrides(self, tmp_path):
        parameters = InputParameters.create("BSFEK")
        parameters.CONTROL.DATASET = "Fe"
        task = TaskDefinition(
            "BSF",
            name="CUSTOM",
            input_files=("{dataset}.aux",),
            cache_files=("{dataset}.cache",),
        )
        registry = InputParameters.input_parameters_registry

        assert registry.files_for_task(
            task, parameters, str(tmp_path), "input"
        ) == (str(tmp_path / "Fe.aux"),)
        assert registry.files_for_task(
            task, parameters, str(tmp_path), "cache"
        ) == (
            str(tmp_path / "Fe.cache"),
            str(tmp_path / "Fe.tau.dir"),
        )

    def test_task_definitions_do_not_own_registries_and_are_pickleable(self):
        task_registry = InputParameters.task_registry
        input_registry = InputParameters.input_parameters_registry
        for name in task_registry.names():
            task = task_registry[name]
            assert not hasattr(task, "_input_parameters_registry")

            restored = pickle.loads(pickle.dumps(task))
            assert restored == task
            assert input_registry.definition_for_task(
                restored
            ) is input_registry.definition_for_task(task)
            restored_parameters = input_registry.create_for_task(restored)
            assert restored_parameters.task_name == task.name
            assert (
                restored_parameters._definition
                is input_registry.definition_for_task(task)
            )

    def test_input_parameters_delegate_file_manifests(self, tmp_path, monkeypatch):
        parameters = InputParameters.create("DOS")
        task = parameters.task_definition
        registry = InputParameters.input_parameters_registry
        calls = []

        def files_for_task(actual_task, actual_parameters, directory, kind):
            calls.append((kind, actual_task, actual_parameters, directory))
            return (f"{kind}.dat",)

        monkeypatch.setattr(registry, "files_for_task", files_for_task)

        assert parameters.input_files(str(tmp_path)) == ("input.dat",)
        assert parameters.output_files(str(tmp_path)) == ("output.dat",)
        assert calls == [
            ("input", task, parameters, str(tmp_path)),
            ("output", task, parameters, str(tmp_path)),
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

    def test_bsfek_kpath_makes_generated_segments_optional_in_input(self):
        ip = InputParameters.create("BSFEK")
        ip.CONTROL.POTFIL = "x"
        ip.TASK.KPATH = 1
        text = ip.to_string(validate=True)
        assert "NKDIR=" not in text
        assert "KA1=" not in text
        assert "KE1=" not in text
        parsed = InputParameters.definition("BSF").read_from_string(text)
        assert parsed.TASK.KPATH() == "1"

    def test_bsfek_explicit_segments_require_input_members(self):
        ip = InputParameters.create("BSFEK")
        ip.CONTROL.POTFIL = "x"
        ip.TASK.NKDIR = 1
        text = ip.to_string(validate=True)
        definition = InputParameters.definition("BSF")
        definition.read_from_string(text)

        def omit(name):
            return "".join(
                line for line in text.splitlines(keepends=True)
                if not line.lstrip().startswith(f"{name}=")
            )

        with pytest.warns(DataValidityError, match="Missing non-optional input item TASK.NKDIR"):
            definition.read_from_string(omit("NKDIR"))
        with pytest.warns(DataValidityError, match="Missing non-optional input item TASK.KE"):
            definition.read_from_string(omit("KE1"))

    def test_bsfkk_required_grid_is_checked_after_parse(self):
        ip = InputParameters.create("BSFKK")
        ip.CONTROL.POTFIL = "x"
        self.set_bsfkk_grid(ip)
        definition = InputParameters.definition("BSF")
        text = ip.to_string(validate=True)
        definition.read_from_string(text)
        without_k1 = "".join(
            line for line in text.splitlines(keepends=True)
            if not line.lstrip().startswith("K1=")
        )
        with pytest.warns(DataValidityError, match="Missing non-optional input item TASK.K1"):
            definition.read_from_string(without_k1)

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
        assert bsfkk.task_name == "BSFKK"
        assert bsfkk.ENERGY.NE().tolist() == [1]
        bsfkk.ENERGY.NE = 2
        assert bsfkk.ENERGY.NE().tolist() == [2]
        assert bsfkk.task_name == "BSFEK"

        bsfek = InputParameters.create("BSFEK")
        assert bsfek.task_name == "BSFEK"
        assert bsfek.ENERGY.NE().tolist() == [200]
        bsfek.ENERGY.NE = 1
        assert bsfek.ENERGY.NE().tolist() == [1]
        assert bsfek.task_name == "BSFKK"

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
            # Both accepted input types select the same KPATH value.
            for value in (path, str(path)):
                ip.TASK.KPATH = value
                assert ip.TASK.KPATH() == str(path)

            # Formatting and parsing depend on the chosen path, not on
            # whether its setter received an int or a numeric string.
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
        for name in InputParameters.task_registry.names():
            ip = InputParameters.for_task(name)
            ip.CONTROL.POTFIL = "xxx"
            definition = ip._definition
            try:
                out = ip.to_string(validate=True)
            except DataValidityError:
                if name == "BSFEK":
                    ip.TASK.KPATH = 1
                    out = ip.to_string(validate=True)
                    with pytest.raises(Exception):
                        ip.TASK.NKDIR = 2
                    ip.TASK.KPATH = None
                    ip.TASK.NKDIR = 2
                    out = ip.to_string(validate=True)
                elif name == "BSFKK":
                    self.set_bsfkk_grid(ip)
                    out = ip.to_string(validate=True)
                elif name == "CHIXAS":
                    # CHIXAS intentionally requires physical type-dependent
                    # NEF/CHISPN/CHIVV arrays, which no generic default supplies.
                    out = ip.to_string(validate=False)
                else:
                    raise
            else:
                if name == "BSFEK":
                    raise Exception("This tasks should not be runnable using the defaults argument")

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
    def test_definitions(self):
        path = os.path.join(os.path.dirname(__file__), "../examples")

        def check(parameters, name):
            assert parameters.name == name
            uname = name.upper()
            assert parameters.__class__ is InputParameters
            assert parameters["CONTROL"]["ADSI"]() == uname
            assert parameters.TASK.TASK() == uname

        for filename_entry in os.listdir(path):
            try:
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
                assert parameters.task_definition is not None
                assert (
                    parameters._definition
                    is parameters.input_parameters_registry.definition_for_task(
                        parameters.task_definition
                    )
                )

                if name == "scf":
                    assert "MODE" not in parameters.to_string()
                if name == "arpes":
                    for task in ["AIPES", "SPLEED", "BAND"]:
                        parameters.TASK.TASK = task

            except Exception as exc:
                raise Exception(
                    f'Parsing of "{filename_entry}" failed with the reason: \n {exc}'
                ).with_traceback(exc.__traceback__)

        definition = InputParameters.definitions["PHAGEN"]
        definition["TASK"]["TASK"].default_value = "scf"
        parameters = definition.read_from_file(os.path.join(path, "phagen.in"))
        assert parameters.TASK.TASK() == "PHAGEN"

        with tempfile.TemporaryFile("w+") as fp:
            parameters.save_to_file(fp)
            fp.seek(0)
            out = fp.read()
            assert "TASKPHAGEN" in re.sub(r"\s", "", out)
        definition["TASK"]["TASK"].default_value = "PHAGEN"
