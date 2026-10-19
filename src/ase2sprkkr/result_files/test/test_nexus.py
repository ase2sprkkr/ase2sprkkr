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
