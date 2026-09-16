import os
import tempfile
import re
import pytest

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..input_parameters import InputParameters
    from ...common.warnings import DataValidityError, DataValidityWarning


class TestDefinitions(TestCase):
    @staticmethod
    def set_bsfkk_grid(ip):
        ip.TASK.NK1 = 10
        ip.TASK.K1 = [1.0, 0.0, 0.0]
        ip.TASK.NK2 = 10
        ip.TASK.K2 = [0.0, 1.0, 0.0]

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
        definition = InputParameters.definition("BSF")

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
        for path in (6, 7, 10):
            ip.TASK.KPATH = path
        with pytest.raises(DataValidityError, match="one of 1-7 or 10"):
            ip.TASK.KPATH = 8

        ip.TASK.KPATH = None
        with pytest.raises(DataValidityError, match="NKDIR cannot be greater than 9"):
            ip.TASK.NKDIR = 10

    def change_task(self):
        ip = InputParameters.create_task("DOS")
        ip.ENERGY.EMIN = 2
        ip.ENERGY.EMAX = 5
        ip.change_task("SCF")
        assert ip.ENERGY.EMIN == 2
        with pytest.raises(AttributeError):
            ip.ENERGY.EMIN

    def jxc(self):
        ip = InputParameters.create_task("SCF")
        ip.MODE.MODE = "nrel"
        with pytest.warns(UserWarning):
            ip.change_task("jxc")
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

                if name == "scf":
                    self.assertFalse("MODE" in t.to_string())
                if name == "arpes":
                    for i in ["AIPES", "SPLEED", "BAND"]:
                        t.TASK.TASK = i

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
