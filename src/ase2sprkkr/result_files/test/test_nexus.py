"""NeXus ARPES exporter regressions (without an SPR-KKR binary)."""
import numpy as np
import pytest
from pathlib import Path

h5py = pytest.importorskip("h5py")
from ase2sprkkr.result_files.nexus import write_arpes_nexus


class SPC:
    def __init__(self, irregular=False):
        self.raw = np.array([[10,-1,100,60,40,.2,1,11],
                             [20,-1,101,62,39,.22,2,12],
                             [10, 0,102,63,39,.24,1,13],
                             [20, 0,103,64,39,.26,2,14]], dtype=float)
        if irregular:
            self.raw[1,1] = -.7

    def RAW_DATA(self): return self.raw
    def NT(self): return 2
    def NP(self): return 1
    def COMMENT(self): return "test"
    def THETA(self): return self.raw[:,0].reshape(2,2)
    def ENERGY(self): return self.raw[:,1].reshape(2,2)
    def TOTAL(self): return self.raw[:,2].reshape(2,2)
    def UP(self): return self.raw[:,3].reshape(2,2)
    def DOWN(self): return self.raw[:,4].reshape(2,2)
    def POLARIZATION(self): return self.raw[:,5].reshape(2,2)
    def K(self): return self.raw[:,6].reshape(2,2)
    def DETERMINANT(self): return self.raw[:,7].reshape(2,2)


def test_export(tmp_path):
    spc = SPC()
    path = tmp_path / "arpes.nxs"
    assert write_arpes_nexus(spc,path) == path
    with h5py.File(path) as f:
        assert f.attrs["NX_class"] == "NXroot"
        assert f["entry"].attrs["NX_class"] == "NXentry"
        assert "definition" not in f["entry"]
        assert f["entry/instrument/electronanalyzer"].attrs["NX_class"] == "NXelectronanalyzer"
        assert "energy_dispersion" not in f["entry/instrument/electronanalyzer"]
        data = f["entry/data"]
        assert list(data.attrs["axes"]) == ["energy", "theta"]
        np.testing.assert_array_equal(data["energy"], [-1,0])
        np.testing.assert_array_equal(data["theta"], [10,20])
        assert data["energy"].attrs["units"] == "eV"
        assert data["theta"].attrs["units"] == "degree"
        assert data["k_parallel"].attrs["units"] == "1/angstrom"
        assert list(data.attrs["k_parallel_indices"]) == [0,1]
        assert "kinetic_energy" not in data
        assert f["entry/instrument/source"].attrs["NX_class"] == "NXsource"
        assert f["entry/instrument/source/probe"].asstr()[()] == "photon"
        assert "energy" not in f["entry/instrument/source"]
        assert f["entry/sample"].attrs["NX_class"] == "NXsample"
        assert "surface_orientation" not in f["entry/sample"]
        assert data["spin_polarization"].attrs["formula"] == "(I_1 - I_2) / (I_1 + I_2)"
        for name, method in (("data",spc.TOTAL),("spin_up",spc.UP),
                             ("spin_down",spc.DOWN),("spin_polarization",spc.POLARIZATION),
                             ("k_parallel",spc.K),("determinant",spc.DETERMINANT)):
            np.testing.assert_array_equal(data[name],method())
        np.testing.assert_array_equal(f["entry/sprkkr_spc/raw_data"],spc.RAW_DATA())


def test_irregular_grid(tmp_path):
    spc = SPC(irregular=True)
    path = tmp_path / "irregular.nxs"
    write_arpes_nexus(spc,path)
    with h5py.File(path) as f:
        d = f["entry/data"]
        assert list(d.attrs["axes"]) == ["energy_index","angle_index"]
        np.testing.assert_array_equal(d["energy"],spc.ENERGY())
        assert list(d.attrs["energy_indices"]) == [0,1]


@pytest.mark.parametrize("name",["DOWN","RAW_DATA"])
def test_invalid_shape(tmp_path,name):
    spc = SPC()
    setattr(spc,name,lambda: np.ones((2,3)))
    path = tmp_path / "bad.nxs"
    with pytest.raises(ValueError,match="(inconsistent|eight columns)"):
        write_arpes_nexus(spc,path)
    assert not path.exists()


def test_input_parameters(tmp_path):
    class Input:
        def to_string(self,*,validate):
            assert validate is False
            return "TASK ARPES\n"
    path = tmp_path / "input.nxs"
    write_arpes_nexus(SPC(),path,input_parameters=Input())
    with h5py.File(path) as f:
        assert f["entry/calculation/parameters/sprkkr_input"].asstr()[()] == "TASK ARPES\n"


def test_spc_api(tmp_path):
    from ase2sprkkr.result_files.definitions.spc import ARPESResultFile
    assert ARPESResultFile.to_nexus(SPC(),tmp_path/"spc.nxs") == tmp_path/"spc.nxs"


def test_result_api(tmp_path):
    from ase2sprkkr.outputs.readers.arpes import ArpesResult
    class Input:
        def to_string(self,*,validate): return "TASK ARPES\n"
    class ResultSPC(SPC):
        def to_nexus(self, filename, *, input_parameters=None, input_file=None,
                     metadata=None):
            return write_arpes_nexus(
                self, filename, input_parameters=input_parameters,
                input_file=input_file, metadata=metadata
            )
    class Result:
        spc = ResultSPC()
        _input_parameters = Input()
        input_parameters_file = None
    path = tmp_path/"result.nxs"
    assert ArpesResult.to_nexus(Result(),path) == path
    with h5py.File(path) as f:
        assert f["entry/calculation/parameters/sprkkr_input"].asstr()[()] == "TASK ARPES\n"


def test_physical_metadata(tmp_path):
    class Value:
        def __init__(self, value):
            self.value = value
        def __call__(self):
            return self.value
    class Section:
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, Value(value))
    class Parameters:
        SPEC_PH = Section(EPHOT=25., THETA=180., PHI=0., POL_P="P")
        ENERGY = Section(EWORK_EV=4.2)
        SPEC_EL = Section(SPOL=4, PSPIN=[0,0,1])
        TASK = Section(MILLER_HKL=[0,0,1])
        def to_string(self, *, validate):
            return "TASK ARPES"
    path = tmp_path / "geometry.nxs"
    write_arpes_nexus(SPC(), path, input_parameters=Parameters())
    with h5py.File(path) as f:
        source = f["entry/instrument/source"]
        assert source["energy"][()] == 25.
        assert source["energy"].attrs["units"] == "eV"
        photon = f["entry/instrument/photon_geometry"]
        assert photon["theta"][()] == 180.
        assert photon["polarization"].asstr()[()] == "P"
        np.testing.assert_array_equal(
            f["entry/sample/surface_orientation/miller_hkl"], [0,0,1])
        spin = f["entry/instrument/spin_projection"]
        np.testing.assert_array_equal(spin["pspin"], [0,0,1])
        np.testing.assert_allclose(
            f["entry/data/kinetic_energy"], SPC().ENERGY() + 20.8)
        assert f["entry/data/kinetic_energy"].attrs["units"] == "eV"


def test_azimuthal_scan_not_misinterpreted(tmp_path):
    spc = SPC()
    spc.NP = lambda: 2
    path = tmp_path / "wrong.nxs"
    with pytest.raises(ValueError, match="NP=1"):
        write_arpes_nexus(spc, path)
    assert not path.exists()

def test_k_scan_source_coordinate_definition(tmp_path):
    class Val:
        def __init__(self, value): self.value = value
        def __call__(self): return self.value
    class Options:
        K1 = Val((1., 0.))
        KA = Val((.5, 0.))
        K2 = Val(None)
    class Parameters:
        SPEC_EL = Options()
        def to_string(self, *, validate): return "TASK ARPES"
    path = tmp_path / "k-scan.nxs"
    write_arpes_nexus(SPC(), path, input_parameters=Parameters())
    with h5py.File(path) as f:
        sampling = f["entry/instrument/electronanalyzer/sampling"]
        assert sampling["mode"].asstr()[()] == "parallel_momentum"
        assert sampling["k1"].attrs["units"] == "1/angstrom"
        assert "distance" in f["entry/data/k_parallel"].attrs["coordinate_definition"]
        np.testing.assert_allclose(sampling["k_x"], [.5, 1.5], atol=1.e-8)
        np.testing.assert_allclose(sampling["k_y"], [0., 0.])


def test_nt_mismatch(tmp_path):
    spc = SPC()
    spc.NT = lambda: 3
    path = tmp_path / "wrong.nt.nxs"
    with pytest.raises(ValueError, match="NT"):
        write_arpes_nexus(spc, path)
    assert not path.exists()


def test_metadata_without_input(tmp_path):
    path = tmp_path / "meta.nxs"
    write_arpes_nexus(SPC(), path, metadata={
        "photon_energy": 21.2,
        "work_function": 4.5,
        "surface_hkl": [1, 1, 1],
        "photon_theta": 180.,
    })
    with h5py.File(path) as f:
        np.testing.assert_allclose(
            f["entry/data/kinetic_energy"], SPC().ENERGY() + 16.7
        )
        assert f["entry/instrument/source/energy"][()] == 21.2
        assert list(f["entry/calculation/export_metadata/keys_overridden_or_supplied"].asstr()[:]) == [
            "photon_energy", "photon_theta", "surface_hkl", "work_function"
        ]


def test_export_metadata_overrides_inputs(tmp_path):
    class Value:
        def __call__(self): return 25.
    class Photon:
        EPHOT = Value()
    class Parameters:
        SPEC_PH = Photon()
        def to_string(self, *, validate): return "TASK ARPES"
    path = tmp_path / "override.nxs"
    write_arpes_nexus(SPC(), path, input_parameters=Parameters(),
                      metadata={"photon_energy": 35.})
    with h5py.File(path) as f:
        assert f["entry/instrument/source/energy"][()] == 35.


def test_unknown_metadata_rejected(tmp_path):
    path = tmp_path / "unknown.nxs"
    with pytest.raises(ValueError, match="Unknown"):
        write_arpes_nexus(SPC(), path, metadata={"not_defined": 123})
    assert not path.exists()


def test_more_native_setup_metadata(tmp_path):
    class V:
        def __init__(self, x): self.x = x
        def __call__(self): return self.x
    class Group:
        def __init__(self, **kw):
            for key, val in kw.items(): setattr(self, key, V(val))
    class Input:
        SPEC_EL = Group(THETA=(0., 60.), PHI=(0., 0.), NT=2, NP=1, TYP=2, SPOL=4)
        SPEC_PH = Group(ALQ=45., DELQ=0., NPOL=1, ICIRC=1, IDREH=0)
        ENERGY = Group(IMV_INI_EV=.05, IMV_FIN_EV=2.)
        SPEC_STR = Group(SURF_BAR=(.25, .25))
        def to_string(self, *, validate): return "TASK ARPES"
    path=tmp_path / "native.nxs"
    write_arpes_nexus(SPC(),path,input_parameters=Input())
    with h5py.File(path) as f:
        sampling=f["entry/instrument/electronanalyzer/sampling"]
        np.testing.assert_array_equal(sampling["theta_range"],[0.,60.])
        assert sampling["typ"][()] == 2
        np.testing.assert_array_equal(f["entry/instrument/spin_projection/pspin"],[0.,0.,1.])
        assert f["entry/instrument/photon_geometry/alq"][()] == 45.
        assert f["entry/calculation/imaginary_potential/final_state"].attrs["units"] == "eV"
        np.testing.assert_array_equal(f["entry/sample/barrier_parameters"],[.25,.25])


def test_real_wse2_arpes_export(tmp_path):
    """The repository's actual SPC sample must survive a NeXus round-trip."""
    from pathlib import Path
    from ase2sprkkr.result_files import ResultFile

    source = Path(__file__).resolve().parents[1] / "examples" / "WSe2_ARPES.spc"
    assert source.is_file(), "Bundled real SPR-KKR ARPES SPC fixture is missing"
    spc = ResultFile.from_file(str(source), try_only="spc")
    energy = np.asarray(spc.ENERGY())
    intensity = np.asarray(spc.TOTAL())
    assert energy.shape == (200, 160)
    target = tmp_path / "WSe2_ARPES.nxs"
    assert spc.to_nexus(target) == target

    with h5py.File(target) as f:
        data = f["entry/data"]
        assert data.attrs["signal"] == "data"
        np.testing.assert_array_equal(data["data"], intensity)
        np.testing.assert_array_equal(f["entry/sprkkr_spc/raw_data"], spc.RAW_DATA())
        assert data["energy"].attrs["units"] == "eV"
        assert data["k_parallel"].attrs["units"] == "1/angstrom"
        assert "definition" not in f["entry"]


def test_independent_nexus_reader(tmp_path):
    """Verify the default plot using the external nexusformat reader."""
    from nexusformat.nexus import nxopen

    target = tmp_path / "external-reader.nxs"
    write_arpes_nexus(SPC(), target)
    with nxopen(str(target)) as root:
        data = root.get_default()
        assert data.nxsignal.shape == (2, 2)
        np.testing.assert_array_equal(data.nxsignal.nxvalue, SPC().TOTAL())


def test_native_photon_polarization_overrides_low_level_options(tmp_path):
    class Opt:
        def __init__(self, value): self.value=value
        def __call__(self): return self.value
    class Section:
        def __init__(self, **items):
            for k, v in items.items(): setattr(self, k, Opt(v))
    class Input:
        SPEC_PH = Section(POL_P="C-", NPOL=1, ICIRC=1, IDREH=0,
                          ALQ=20., DELQ=0., THETA=45., EPHOT=25.)
        def to_string(self, *, validate): return "POL_P=C-"
    path = tmp_path/"c_minus.nxs"
    write_arpes_nexus(SPC(), path, input_parameters=Input())
    with h5py.File(path) as f:
        g=f["entry/instrument/photon_geometry"]
        assert g["npol"][()] == 2
        assert g["icirc"][()] == 2
        assert g["idreh"][()] == -1
        assert g["alq"][()] == 45.
        assert g["delq"][()] == 90.


def test_native_variable_photon_incidence(tmp_path):
    path=tmp_path/"variable.nxs"
    write_arpes_nexus(SPC(),path,metadata={
        "photon_polarization":"P", "photon_ifsp":1,
        "photon_theta_fix":15., "photon_theta":45.
    })
    with h5py.File(path) as f:
        g=f["entry/instrument/photon_geometry"]
        assert g["ifsp"][()] == 1
        assert g["theta_fix"][()] == 15.
        assert "theta" not in g
        np.testing.assert_array_equal(g["polar_angle_before_tilt"],SPC().THETA()-15.)


def test_external_legacy_input_does_not_claim_modern_geometry(tmp_path):
    class Opt:
        def __init__(self,v):self.v=v
        def __call__(self):return self.v
    class Section:
        def __init__(self,**options):
            for key,val in options.items():setattr(self,key,Opt(val))
    class Input:
        TASK=Section(INPVER=0,STRVER=0,MILLER_HKL=[0,0,1])
        SPEC_PH=Section(POL_P="P",THETA=30.,EPHOT=25.)
        ENERGY=Section(EWORK_EV=4.2)
        def to_string(self,*,validate):return "INPVER=0"
    path=tmp_path/"legacy.nxs"
    write_arpes_nexus(SPC(),path,input_parameters=Input())
    with h5py.File(path) as f:
        assert "energy" not in f["entry/instrument/source"]
        assert "surface_orientation" not in f["entry/sample"]
        assert "kinetic_energy" not in f["entry/data"]
        g=f["entry/calculation/external_native_input"]
        assert not bool(g["geometry_complete"][()])
        assert g["input_version"][()] == 0


def test_native_spin_polarization_matrix(tmp_path):
    class Opt:
        def __init__(self,v):self.v=v
        def __call__(self):return self.v
    class Group:
        def __init__(self,**kwargs):
            for name,v in kwargs.items():setattr(self,name,Opt(v))
    class Input:
        SPEC_EL=Group(SPOL=4,POL_E="PX",
                      Q1=1.+0.j,Q2=0.j,Q3=0.j,Q4=1.+0.j,
                      POL0=[1,0,0],POL0L=[0,1,0])
        def to_string(self,*,validate):return "POL_E=PX"
    path=tmp_path/"spin.nxs"
    write_arpes_nexus(SPC(),path,input_parameters=Input())
    with h5py.File(path) as f:
        spin=f["entry/instrument/spin_projection"]
        np.testing.assert_allclose(
            spin["q_amplitudes"], [.707,.707,.707,-.707]
        )
        assert "pol0" not in spin
        assert "pol0l" not in spin



def test_warns_when_input_unavailable(tmp_path):
    target = tmp_path / "no-input.nxs"
    with pytest.warns(UserWarning, match="No SPR-KKR input parameters"):
        write_arpes_nexus(SPC(), target)
    with h5py.File(target) as f:
        assert "parameters" not in f["entry/calculation"]
        assert "energy" not in f["entry/instrument/source"]


def test_explicit_input_file_wins_over_parameters(tmp_path, monkeypatch):
    from ase2sprkkr.input_parameters.input_parameters import InputParameters

    source = tmp_path / "explicit.inp"
    source.write_text("dummy", encoding="utf-8")
    class Opt:
        def __call__(self): return 48.
    class Group:
        EPHOT = Opt()
    class Parsed:
        SPEC_PH = Group()
        def to_string(self, *, validate): return "PARSED INPUT"
    parsed = Parsed()
    calls = []
    def read(cls, filename):
        calls.append(Path(filename))
        return parsed
    monkeypatch.setattr(InputParameters, "from_file", classmethod(read))

    target = tmp_path / "explicit.nxs"
    write_arpes_nexus(SPC(), target, input_file=source,
                      input_parameters=object())
    assert calls == [source]
    with h5py.File(target) as f:
        assert f["entry/instrument/source/energy"][()] == 48.
        assert f["entry/calculation/parameters/source_file"].asstr()[()] == str(source)


def test_companion_autodiscovery(tmp_path, monkeypatch):
    from pathlib import Path
    from ase2sprkkr.input_parameters.input_parameters import InputParameters

    source = tmp_path / "WSe2_ARPES_data.spc"
    source.touch()
    companion = tmp_path / "WSe2.inp"
    companion.write_text("dummy", encoding="utf-8")
    spc = SPC()
    spc._filename = str(source)

    class Parsed:
        def to_string(self, *, validate): return "FOUND INPUT"
    monkeypatch.setattr(InputParameters, "from_file",
                        classmethod(lambda cls, file: Parsed()))
    target = tmp_path / "auto.nxs"
    write_arpes_nexus(spc, target)
    with h5py.File(target) as f:
        assert f["entry/calculation/parameters/sprkkr_input"].asstr()[()] == "FOUND INPUT"
        assert f["entry/calculation/parameters/source_file"].asstr()[()] == str(companion)


def test_ambiguous_companion_does_not_pick_arbitrarily(tmp_path):
    (tmp_path / "first.inp").touch()
    (tmp_path / "second.inp").touch()
    spc = SPC()
    spc._filename = str(tmp_path / "some_ARPES_data.spc")
    with pytest.warns(UserWarning, match="No SPR-KKR input parameters"):
        write_arpes_nexus(spc, tmp_path / "uncertain.nxs")


def test_explicit_missing_input_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="input file not found"):
        write_arpes_nexus(
            SPC(), tmp_path / "out.nxs", input_file=tmp_path / "missing.inp"
        )
    assert not (tmp_path / "out.nxs").exists()


def test_output_file_retains_source_filename():
    from pathlib import Path
    from ase2sprkkr.result_files import ResultFile
    path = Path(__file__).resolve().parents[1] / "examples" / "WSe2_ARPES.spc"
    spc = ResultFile.from_file(str(path), try_only="spc")
    assert Path(spc._filename) == path


def test_result_passes_explicit_input_and_input_parameters():
    from ase2sprkkr.outputs.readers.arpes import ArpesResult

    class StubSPC:
        def to_nexus(self, target, **kwargs):
            return target, kwargs

    class Dummy:
        _input_parameters = "already-known"
        directory = "/tmp/calc"
        input_parameters_file = None
        spc = StubSPC()

    target, kwargs = ArpesResult.to_nexus(
        Dummy(), "output.nxs", input_file="other.inp"
    )
    assert target == "output.nxs"
    assert Path(kwargs["input_file"]) == Path("/tmp/calc/other.inp")
    assert kwargs["input_parameters"] == "already-known"


def test_kxky_spc_mode_is_not_mistaken_for_energy_spectrum(tmp_path):
    """The current native SPC parser distinguishes energy and kx/ky modes."""
    class KxKy(SPC):
        def MODE(self):
            return "kx_ky"
    target = tmp_path / "kxky.nxs"
    with pytest.raises(ValueError, match="kx_ky"):
        write_arpes_nexus(KxKy(), target)
    assert not target.exists()
