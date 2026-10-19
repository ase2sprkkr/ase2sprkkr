"""Regression tests of native SPR-KKR SPC → NeXus conversion."""
from pathlib import Path
import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
from ase2sprkkr.result_files.nexus import write_arpes_nexus


class SPC:
    def __init__(self, mode="energy", np_=1):
        self.mode, self.np_ = mode, np_
        self.nt, self.ne = 3, 2
        rows = []
        energies = range(self.ne) if mode == "energy" else range(1)
        for ie in energies:
            for jp in range(np_):
                for jt in range(self.nt):
                    if mode == "energy":
                        rows.append([10 + jt, -1. + ie, 1 + jt + 3*jp + 6*ie,
                                     1., 2., -.333, .5+jt, 3.])
                    else:
                        rows.append([.5 + jt + .25*jp, -1. + .5*jp,
                                     1+jt+3*jp, 1.,2.,-.333,3.])
        self.raw=np.asarray(rows,float)

    def MODE(self): return self.mode
    def RAW_DATA(self): return self.raw
    def NT(self): return self.nt
    def NP(self): return self.np_
    def NE(self): return self.ne if self.mode=="energy" else 1
    def COMMENT(self): return "synthetic native spectrum"
    def spc_columns(self):
        from ase2sprkkr.result_files.definitions.spc import SPC_COLUMNS
        return SPC_COLUMNS[self.MODE()]
    def _field(self, col):
        shape = (self.ne, self.np_, self.nt) if self.mode == "energy" and self.np_ > 1 else (
            (self.ne, self.nt) if self.mode == "energy" else (self.np_, self.nt)
        )
        return self.raw[:, col].reshape(shape)
    def TOTAL(self): return self._field(2)
    def UP(self): return self._field(3)
    def DOWN(self): return self._field(4)
    def POLARIZATION(self): return self._field(5)
    def DETERMINANT(self): return self._field(7 if self.mode == "energy" else 6)
    def ENERGY(self): return self._field(1)
    def THETA(self): return self._field(0)
    def K(self): return self._field(6)
    def KX(self): return self._field(0)
    def KY(self): return self._field(1)


@pytest.mark.parametrize("mode,np_", [("energy",1),("energy",2),("kx_ky",2)])
def test_native_modes(tmp_path, mode, np_):
    spc=SPC(mode,np_)
    dest=tmp_path/"map.nxs"
    with pytest.warns(UserWarning,match="No SPR-KKR input"):
        assert write_arpes_nexus(spc,dest)==dest
    expected=(2,2,3) if mode=="energy" and np_>1 else (
        (2,3) if mode=="energy" else (np_,3)
    )
    with h5py.File(dest) as f:
        data=f["entry/data"]
        assert data["data"].shape==expected
        np.testing.assert_array_equal(data["data"],spc.raw[:,2].reshape(expected))
        np.testing.assert_array_equal(f["entry/sprkkr_spc/raw_data"],spc.raw)
        assert len(f["entry/sprkkr_spc/raw_data"].attrs["columns"])==spc.raw.shape[1]
        assert "original_input" not in f["entry/calculation"]
        if mode=="kx_ky":
            assert "energy" not in data
            assert "coordinate_warning" in data.attrs
            assert list(data.attrs["axes"]) == ["k_y", "k_x"]
            np.testing.assert_allclose(data["k_x"], spc.KX()[0])
            np.testing.assert_allclose(data["k_y"], spc.KY()[:, 0])
            assert "k_origin" not in data
            assert "k_basis" not in data
        else:
            assert data["energy"].attrs["units"]=="eV"


def test_explicit_input_precedence_and_original_options(tmp_path,monkeypatch):
    from ase2sprkkr.input_parameters.input_parameters import InputParameters
    source=tmp_path/"correct.inp"
    source.write_text("TASK ARPES\n",encoding="utf-8")
    class Opt:
        def __init__(self,value): self.value=value
        def __call__(self): return self.value
    class Sect:
        def __init__(self,**kw):
            for k,v in kw.items():setattr(self,k,Opt(v))
    class Parameters:
        SPEC_PH=Sect(EPHOT=25.,POL_P="P",ALQ=45.)
        SPEC_EL=Sect(KA=[.2,.3],K1=[1,0],K2=[0,1],NK1=3,NK2=2)
        TASK=Sect(INPVER=1,STRVER=1,MILLER_HKL=[0,0,1])
        def to_string(self,*,validate): return "WRONG OBJECT"
    def read(cls,p):
        assert Path(p)==source
        return Parameters()
    monkeypatch.setattr(InputParameters,"from_file",classmethod(read))
    out=tmp_path/"case.nxs"
    write_arpes_nexus(SPC("kx_ky",2),out,
                      input_file=source,input_parameters=object())
    with h5py.File(out) as f:
        parameters=f["entry/calculation/parameters"]
        assert parameters["sprkkr_input"].asstr()[()]=="TASK ARPES\n"
        original=f["entry/calculation/original_input"]
        np.testing.assert_array_equal(original["SPEC_EL/NK1"],3)
        np.testing.assert_array_equal(original["SPEC_EL/K1"],[1,0])
        assert "k_basis" not in f["entry/data"]


def test_missing_input_warning(tmp_path):
    with pytest.warns(UserWarning,match="No SPR-KKR input"):
        write_arpes_nexus(SPC(),tmp_path/"missing.nxs")


def test_explicit_invalid_input_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        write_arpes_nexus(SPC(),tmp_path/"no.nxs",
                          input_file=tmp_path/"missing.inp")
    assert not (tmp_path/"no.nxs").exists()


def test_companion_autodetection(tmp_path,monkeypatch):
    from ase2sprkkr.input_parameters.input_parameters import InputParameters
    p=tmp_path/"Run_ARPES_data.spc"
    p.touch()
    inp=tmp_path/"Run.inp"
    inp.write_text("TASK ARPES\n")
    spc=SPC()
    spc._filename=str(p)
    class Parameters:
        def to_string(self,*,validate): return "TASK ARPES"
    monkeypatch.setattr(InputParameters,"from_file",classmethod(
        lambda cls,f: Parameters()))
    out=tmp_path/"detected.nxs"
    write_arpes_nexus(spc,out)
    with h5py.File(out) as f:
        assert f["entry/calculation/parameters/source_file"].asstr()[()]==str(inp)


def test_ambiguous_input_not_guessed(tmp_path):
    p=tmp_path/"Run_ARPES_data.spc"
    p.touch()
    for suffix in [".in",".inp"]:
        (tmp_path/("Run"+suffix)).touch()
    spc=SPC()
    spc._filename=str(p)
    with pytest.warns(UserWarning,match="No SPR-KKR input"):
        write_arpes_nexus(spc,tmp_path/"out.nxs")


def test_atomic_output_on_failure(tmp_path):
    target=tmp_path/"existing.nxs"
    target.write_bytes(b"do not overwrite")
    invalid=SPC()
    invalid.raw=np.ones((4,7))
    with pytest.raises(ValueError,match="8 columns"):
        write_arpes_nexus(invalid,target)
    assert target.read_bytes()==b"do not overwrite"
    assert list(tmp_path.iterdir())==[target]


def test_real_wse2_spectrum(tmp_path):
    from ase2sprkkr.result_files import ResultFile
    source=Path(__file__).resolve().parents[1]/"examples"/"WSe2_ARPES.spc"
    spc=ResultFile.from_file(str(source),try_only="spc")
    out=tmp_path/"real.nxs"
    with pytest.warns(UserWarning,match="No SPR-KKR input"):
        spc.to_nexus(out)
    with h5py.File(out) as f:
        np.testing.assert_array_equal(f["entry/data/data"],spc.TOTAL())
        np.testing.assert_array_equal(f["entry/sprkkr_spc/raw_data"],spc.RAW_DATA())


def test_independent_nexus_reader(tmp_path):
    from nexusformat.nexus import nxopen
    dest=tmp_path/"reader.nxs"
    with pytest.warns(UserWarning):
        write_arpes_nexus(SPC("kx_ky",2),dest)
    with nxopen(str(dest)) as root:
        assert root.get_default().nxsignal.shape==(2,3)



def test_atomic_output_filename_is_final_target(tmp_path):
    target=tmp_path/"real-target.nxs"
    with pytest.warns(UserWarning):
        write_arpes_nexus(SPC(),target)
    with h5py.File(target) as f:
        assert f.attrs["file_name"]=="real-target.nxs"


def test_native_au_constant_energy_spc(tmp_path):
    """The repository's actual seven-column kx/ky calculation must export."""
    from ase2sprkkr.result_files import ResultFile
    path=Path(__file__).resolve().parents[1]/"examples"/"Au_kx_ky_ARPES_data.spc"
    spc=ResultFile.from_file(str(path),try_only="spc")
    assert spc.MODE()=="kx_ky"
    assert spc.RAW_DATA().shape[1]==7
    with pytest.warns(UserWarning,match="No SPR-KKR input"):
        spc.to_nexus(tmp_path/"au_kmap.nxs")
    with h5py.File(tmp_path/"au_kmap.nxs") as f:
        assert f["entry/data"].attrs["sprkkr_mode"]=="kx_ky"
        np.testing.assert_array_equal(f["entry/data/data"],spc.TOTAL())
        np.testing.assert_array_equal(f["entry/sprkkr_spc/raw_data"],spc.RAW_DATA())


def test_no_physics_reconstructed_from_input(tmp_path):
    """Original vectors are provenance; export must not recalculate K coordinates."""
    class Value:
        def __init__(self,value): self.value=value
        def __call__(self): return self.value
    class Group:
        def __init__(self,**values):
            for name,value in values.items(): setattr(self,name,Value(value))
    class Parameters:
        SPEC_EL=Group(KA=[.3,.4],K1=[1.5,0.],K2=[.5,2.],NK1=3,NK2=2)
        def to_string(self,*,validate): return "TASK ARPES"
    path=tmp_path/"native.nxs"
    write_arpes_nexus(SPC("kx_ky",2),path,input_parameters=Parameters())
    with h5py.File(path) as f:
        spec=f["entry/calculation/original_input/SPEC_EL"]
        np.testing.assert_array_equal(spec["KA"],[.3,.4])
        np.testing.assert_array_equal(spec["K1"],[1.5,0.])
        np.testing.assert_array_equal(spec["K2"],[.5,2.])
        assert spec["NK1"][()]==3
        assert spec["NK2"][()]==2
        data=f["entry/data"]
        assert "k_basis" not in data
        assert data.attrs["scan_definition"]=="../calculation/original_input/SPEC_EL"



def test_arpes_static_azimuth_selection(tmp_path, monkeypatch):
    """Multiple file slices and live mode must share the same plotting path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ase2sprkkr.result_files.definitions.spc import ARPESResultFile
    from ase2sprkkr.result_files import ResultFile
    from ase2sprkkr.gui.plot import colormesh

    source = Path(__file__).resolve().parents[1] / "examples" / "WSe2_ARPES.spc"
    # Use a minimal adapted object to exercise the NP>1 branch without
    # creating a large file or launching a GUI.
    class Value:
        def __init__(self, value, name): self.value, self.name = value, name
        def __call__(self): return self.value
        def plot(self, **kwargs): return kwargs
    class PlotResult:
        MODE = lambda self: "energy"
        NP = lambda self: 5
        TOTAL = Value(np.zeros((2, 5, 3)), "TOTAL")
        UP = Value(np.zeros((2, 5, 3)), "UP")
        DOWN = Value(np.zeros((2, 5, 3)), "DOWN")
        POLARIZATION = Value(np.zeros((2, 5, 3)), "POLARIZATION")

    calls = []
    def record(option, **kwargs):
        calls.append(kwargs)
    for val in (PlotResult.TOTAL, PlotResult.UP, PlotResult.DOWN, PlotResult.POLARIZATION):
        monkeypatch.setattr(val, "plot", lambda **kwargs: calls.append(kwargs))
    ARPESResultFile.plot(PlotResult(), filename=tmp_path / "one.png", live=False)
    assert len(calls) == 4
    assert all(item["azimuth_indices"] == (0,) for item in calls)
    calls.clear()
    ARPESResultFile.plot(
        PlotResult(), filename=tmp_path / "three.png", live=False, azimuths=3
    )
    assert len(calls) == 12
    assert [calls[i]["azimuth_indices"] for i in (0, 4, 8)] == [(0,), (2,), (4,)]
    calls.clear()
    ARPESResultFile.plot(
        PlotResult(), filename=tmp_path / "angles.png",
        live=False, azimuths=[15., 45.], phi_range=(0., 60.)
    )
    assert [calls[i]["azimuth_indices"] for i in (0, 4)] == [(1,), (3,)]
    calls.clear()
    ARPESResultFile.plot(PlotResult(), show=False, live=True)
    assert len(calls) == 4
    assert all(item["azimuth_indices"] == (0, 1, 2, 3, 4) for item in calls)
    plt.close("all")


def test_colormesh_callable_slider_shared():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ase2sprkkr.gui.plot import colormesh
    fig, axes = plt.subplots(1, 2)
    xx, yy = np.meshgrid(np.arange(3), np.arange(3))
    calls = []
    for ax in axes:
        colormesh(xx, yy, lambda i: calls.append(i) or np.full((3, 3), i + 1.),
                  slider_values=(0, 1, 2), live=True, axis=ax)
    widget = fig._ase2sprkkr_slider[0]
    assert calls == [0, 0]
    widget.set_val(2)
    assert calls == [0, 0, 2, 2]
    assert len(fig.axes) == 3  # one common slider, not one per subplot
    plt.close(fig)
