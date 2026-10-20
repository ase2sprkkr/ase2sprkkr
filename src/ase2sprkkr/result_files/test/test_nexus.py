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
                        rows.append([.5 + jt, -1. + .5*jp,
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
            assert "coordinate_warning" not in data.attrs
            assert list(data.attrs["axes"]) == ["k_y", "k_x"]
            np.testing.assert_allclose(data["k_x"], spc.KX()[0])
            np.testing.assert_allclose(data["k_y"], spc.KY()[:, 0])
            assert data["k_x"].attrs["units"] == "1/angstrom"
            assert data["k_y"].attrs["units"] == "1/angstrom"
            assert "k_origin" not in data
            assert "k_basis" not in data
        else:
            assert data["energy"].attrs["units"]=="eV"
            assert data["k_parallel"].attrs["units"] == "1/angstrom"
            assert "native_units" not in data["k_parallel"].attrs


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
    nxopen = pytest.importorskip("nexusformat.nexus").nxopen
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
    # Use a minimal adapted object to exercise the NP>1 branch without
    # creating a large file or launching a GUI.
    class Value:
        def __init__(self, value, name): self.value, self.name = value, name
        def __call__(self): return self.value
        def plot(self, **kwargs): return kwargs
    class PlotResult:
        MODE = lambda self: "energy"
        NP = lambda self: 5
        NT = lambda self: 3
        TOTAL = Value(np.zeros((2, 5, 3)), "TOTAL")
        UP = Value(np.zeros((2, 5, 3)), "UP")
        DOWN = Value(np.zeros((2, 5, 3)), "DOWN")
        POLARIZATION = Value(np.zeros((2, 5, 3)), "POLARIZATION")

    calls = []
    def record(option, **kwargs):
        calls.append(kwargs)
    for val in (PlotResult.TOTAL, PlotResult.UP, PlotResult.DOWN, PlotResult.POLARIZATION):
        monkeypatch.setattr(val, "plot", lambda **kwargs: calls.append(kwargs))
    ARPESResultFile.plot(PlotResult(), filename=tmp_path / "one.png", dpi=50, live=False)
    assert len(calls) == 4
    assert all(item["azimuth_indices"] == (0,) for item in calls)
    calls.clear()
    ARPESResultFile.plot(
        PlotResult(), filename=tmp_path / "three.png", dpi=50, live=False, azimuths=3
    )
    assert len(calls) == 12
    assert [calls[i]["azimuth_indices"] for i in (0, 4, 8)] == [(0,), (2,), (4,)]
    calls.clear()
    ARPESResultFile.plot(
        PlotResult(), filename=tmp_path / "angles.png", dpi=50,
        live=False, azimuths=[15., 45.], phi_range=(0., 60.)
    )
    assert [calls[i]["azimuth_indices"] for i in (0, 4)] == [(1,), (3,)]
    calls.clear()
    ARPESResultFile.plot(PlotResult(), show=False, live=True)
    assert len(calls) == 4
    assert all(item["azimuth_indices"] == (0, 1, 2, 3, 4) for item in calls)
    plt.close("all")


def test_colormesh_callable_controls_shared():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ase2sprkkr.gui.plot import colormesh, PlotControls, Multiplot
    xx, yy = np.meshgrid(np.arange(3), np.arange(3))
    calls = []
    with Multiplot(layout=(1, 2), show=False, live=True,
                   live_options={"direction": {"options": ("NP", "NT")},
                                 "slice": {"slider": lambda state: range(
                                     3 if state["direction"] == "NP" else 5)}}) as mp:
        for _ in range(2):
            with mp.new_axis() as (axis, kwargs):
                colormesh(xx, yy, lambda state: calls.append(dict(state)) or
                          np.full((3, 3), state["slice"] + 1),
                          axis=axis, **kwargs)
    controls = mp.controls
    assert calls == [{"direction": "NP", "slice": 0}] * 2
    controls.widgets["slice"].set_val(2)
    assert calls[-1] == {"direction": "NP", "slice": 2}
    controls.widgets["direction"].set_active(1)
    assert controls.state == {"direction": "NT", "slice": 2}
    assert controls._slider_values["slice"] == tuple(range(5))
    controls.widgets["slice"].set_val(4)
    assert calls[-1] == {"direction": "NT", "slice": 4}
    assert len(mp.figure.axes) == 4  # 2 panels, options and slider
    plt.close(mp.figure)


def test_standalone_controls_and_static():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ase2sprkkr.gui.plot import colormesh
    xx, yy = np.meshgrid(np.arange(3), np.arange(3))
    opts = {"slice": {"slider": range(3)}}
    colormesh(xx, yy, lambda state: np.full((3, 3), state["slice"] + 1),
              live=True, live_options=opts, show=False)
    fig = plt.gcf()
    assert fig._ase2sprkkr_controls.widgets["slice"] is not None
    plt.close(fig)
    colormesh(xx, yy, np.ones((3, 3)), show=False)
    assert not hasattr(plt.gcf(), "_ase2sprkkr_controls")
    plt.close("all")


def test_static_separate_plots_with_callable_data(tmp_path):
    """Shared initial live_options state works without a visible figure/widget."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from ase2sprkkr.gui.plot import Multiplot, colormesh, lineplot
    xx, yy = np.meshgrid(np.arange(3), np.arange(3))
    states = []
    filename = tmp_path / "section.png"
    with Multiplot(
        layout=(1, 2), separate_plots=True, live=False,
        live_options={"slice": {"slider": (1, 2, 3)}},
        filename=filename, show=False, dpi=50,
    ) as mp:
        with mp.new_axis(name="map") as (ax, kwargs):
            colormesh(
                xx, yy,
                lambda state: states.append(state.copy()) or
                np.full((3, 3), state["slice"]),
                axis=ax, **kwargs,
            )
        with mp.new_axis(name="line") as (ax, kwargs):
            lineplot([0, 1], lambda state: [state["slice"], 2], axis=ax, **kwargs)
    assert states == [{"slice": 1}]
    assert (tmp_path / "section_map.png").is_file()
    assert (tmp_path / "section_line.png").is_file()
    assert mp.controls.state == {"slice": 1}
    assert not mp.controls.widgets
    plt.close("all")


def test_nexus_axes_are_in_storage_order_not_plot_xy_order(tmp_path):
    """For a 2D image, signal[0, :] is the first y row, not the first x column."""
    spc = SPC("kx_ky", 2)
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        write_arpes_nexus(spc, tmp_path / "orientation.nxs")
    with h5py.File(tmp_path / "orientation.nxs") as f:
        group = f["entry/data"]
        assert list(group.attrs["axes"]) == ["k_y", "k_x"]
        assert group.attrs["k_y_indices"] == 0
        assert group.attrs["k_x_indices"] == 1
        assert group["data"][1, 2] == spc.raw[1 * spc.NT() + 2, 2]
        assert group["k_x"][2] == spc.raw[2, 0]
        assert group["k_y"][1] == spc.raw[spc.NT(), 1]



def test_energy_theta_axes_follow_storage_order(tmp_path):
    # Intensity has rows of constant energy, with theta varying by column.
    spc = SPC("energy")
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        write_arpes_nexus(spc, tmp_path / "energy_axes.nxs")
    with h5py.File(tmp_path / "energy_axes.nxs") as f:
        group = f["entry/data"]
        assert list(group.attrs["axes"]) == ["energy", "theta"]
        assert group.attrs["energy_indices"] == 0
        assert group.attrs["theta_indices"] == 1
        np.testing.assert_array_equal(group["energy"], spc.ENERGY()[:, 0])
        np.testing.assert_array_equal(group["theta"], spc.THETA()[0, :])
        assert group["data"][1, 2] == spc.raw[spc.NT() + 2, 2]


def test_oblique_k_grid_uses_native_two_dimensional_coordinates(tmp_path):
    spc = SPC("kx_ky", 3)
    for i in range(spc.NP()):
        for j in range(spc.NT()):
            n = i * spc.NT() + j
            spc.raw[n, 0] = j + 0.3 * i
            spc.raw[n, 1] = i + 0.2 * j
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        write_arpes_nexus(spc, tmp_path / "oblique.nxs")
    with h5py.File(tmp_path / "oblique.nxs") as f:
        group = f["entry/data"]
        assert list(group.attrs["axes"]) == ["k_y", "k_x"]
        assert list(group.attrs["k_x_indices"]) == [0, 1]
        assert list(group.attrs["k_y_indices"]) == [0, 1]
        assert "coordinate_warning" not in group.attrs
        assert group["k_x"].attrs["units"] == "1/angstrom"
        assert group["k_y"].attrs["units"] == "1/angstrom"
        np.testing.assert_array_equal(group["k_x"], spc.KX())
        np.testing.assert_array_equal(group["k_y"], spc.KY())
        assert group["k_x"][1, 1] == pytest.approx(1.3)
        assert group["k_y"][1, 1] == pytest.approx(1.2)
        assert group["data"][1, 1] == spc.TOTAL()[1, 1]



def test_small_nonseparable_offset_is_not_lost(tmp_path):
    # The former np.allclose default tolerance would drop this native value.
    spc = SPC("kx_ky", 2)
    spc.raw[spc.NT(), 0] += 1e-10
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        write_arpes_nexus(spc, tmp_path / "small_offset.nxs")
    with h5py.File(tmp_path / "small_offset.nxs") as f:
        group = f["entry/data"]
        assert list(group.attrs["axes"]) == ["k_y", "k_x"]
        # Only k_x is coupled; k_y remains a genuine one-dimensional axis.
        assert group["k_x"].shape == (spc.NP(), spc.NT())
        assert list(group.attrs["k_x_indices"]) == [0, 1]
        assert group["k_y"].shape == (spc.NP(),)
        assert group.attrs["k_y_indices"] == 0
        np.testing.assert_array_equal(group["k_x"], spc.KX())
        assert group["k_x"][1, 0] != group["k_x"][0, 0]


def test_three_dimensional_physical_coordinates_are_default(tmp_path):
    spc = SPC("energy", 2)
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        write_arpes_nexus(spc, tmp_path / "angle3d.nxs")
    with h5py.File(tmp_path / "angle3d.nxs") as f:
        data = f["entry/data"]
        assert list(data.attrs["axes"]) == ["energy", "azimuth_index", "theta"]
        assert data.attrs["energy_indices"] == 0
        assert data.attrs["theta_indices"] == 2
        assert data["energy"].shape == (spc.ne,)
        assert data["theta"].shape == (spc.nt,)
        np.testing.assert_array_equal(data["data"], spc.TOTAL())

    # A coupled coordinate must retain all its physical values instead
    # of being silently replaced by a coordinate index.
    spc.raw[:, 0] += np.repeat(np.arange(spc.ne), spc.np_ * spc.nt)
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        write_arpes_nexus(spc, tmp_path / "coupled3d.nxs")
    with h5py.File(tmp_path / "coupled3d.nxs") as f:
        data = f["entry/data"]
        assert data["theta"].shape == (spc.ne, spc.np_, spc.nt)
        assert list(data.attrs["theta_indices"]) == [0, 1, 2]
        np.testing.assert_array_equal(data["theta"], spc.THETA())


def test_native_k_parallel_is_inverse_angstrom(tmp_path):
    from ase2sprkkr.result_files import ResultFile
    source = Path(__file__).resolve().parents[1] / "examples" / "WSe2_ARPES.spc"
    spc = ResultFile.from_file(str(source), try_only="spc")
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        spc.to_nexus(tmp_path / "units.nxs")
    with h5py.File(tmp_path / "units.nxs") as f:
        ds = f["entry/data/k_parallel"]
        assert ds.attrs["units"] == "1/angstrom"
        assert "native_units" not in ds.attrs
        np.testing.assert_array_equal(ds, spc.K())


def test_extended_scan_provenance(tmp_path):
    class Opt:
        def __init__(self, value): self.value = value
        def __call__(self): return self.value
    class Sect:
        def __init__(self, **options):
            for key, value in options.items():
                setattr(self, key, Opt(value))
    class Parameters:
        TASK = Sect(INPVER=1, STRVER=1, IQ_AT_SURF=2)
        SPEC_EL = Sect(K3=[1., 0.], K4=[0., 1.], NK3=7, NK4=8,
                       TYP=3, ISTR=[1, 2], POL0=[0, 0, 1])
        SPEC_STR = Sect(NLAT_G_VEC=57, TRANSP_BAR=False)
        def to_string(self, *, validate): return "TASK ARPES"
    write_arpes_nexus(SPC(), tmp_path / "provenance.nxs",
                      input_parameters=Parameters())
    with h5py.File(tmp_path / "provenance.nxs") as f:
        original = f["entry/calculation/original_input"]
        assert original["TASK/IQ_AT_SURF"][()] == 2
        assert original["SPEC_EL/NK3"][()] == 7
        assert original["SPEC_EL/NK4"][()] == 8
        assert original["SPEC_EL/TYP"][()] == 3
        np.testing.assert_array_equal(original["SPEC_EL/K3"], [1., 0.])
        np.testing.assert_array_equal(original["SPEC_EL/K4"], [0., 1.])
        assert original["SPEC_STR/NLAT_G_VEC"][()] == 57
        assert not bool(original["SPEC_STR/TRANSP_BAR"][()])


def test_native_wse2_k_parallel_physical_units_from_energy_dispersion():
    """Energy spacing independently establishes 1/angstrom, not pi/angstrom.

    At fixed theta, the free-electron relation yields
    dE = (hbar^2 / 2m) * d(k_parallel^2) / sin(theta)^2.
    This does not require knowing the photon energy or work function.
    """
    from ase2sprkkr.result_files import ResultFile
    source = Path(__file__).resolve().parents[1] / "examples" / "WSe2_ARPES.spc"
    spc = ResultFile.from_file(str(source), try_only="spc")
    assert spc.NP() == 1
    # The NT points at each energy have the same polar angles.
    last_angle_per_energy = spc.RAW_DATA()[spc.NT() - 1::spc.NT()]
    assert np.allclose(last_angle_per_energy[:, 0], 20)
    angle = np.deg2rad(last_angle_per_energy[:, 0])
    k_parallel = last_angle_per_energy[:, 6]
    kinetic_energy = 3.80998212 * (k_parallel / np.sin(angle))**2
    initial_energy = last_angle_per_energy[:, 1]
    slope = np.polyfit(initial_energy, kinetic_energy, 1)[0]
    assert slope == pytest.approx(1.0, abs=0.002)
    # Misinterpreting this momentum as pi/angstrom would give slope pi**2.
    assert slope * np.pi**2 != pytest.approx(1.0, abs=0.002)


def test_native_mo2se4_k_parallel_physical_unit(tmp_path):
    """The native SPC momentum data are exported in inverse angstroms."""
    from ase2sprkkr.result_files import ResultFile
    source = Path(__file__).resolve().parents[1] / "examples" / "Mo2Se4_ARPES_ARPES_data.spc"
    spc = ResultFile.from_file(str(source), try_only="spc")
    with pytest.warns(UserWarning, match="No SPR-KKR input"):
        spc.to_nexus(tmp_path / "mo2se4.nxs")
    with h5py.File(tmp_path / "mo2se4.nxs") as f:
        ds = f["entry/data/k_parallel"]
        assert ds.attrs["units"] == "1/angstrom"
        assert "native_units" not in ds.attrs
        np.testing.assert_array_equal(ds, spc.K())
