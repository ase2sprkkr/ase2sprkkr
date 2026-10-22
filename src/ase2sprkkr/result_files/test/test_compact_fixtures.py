"""Compact real-data-derived smoke fixtures, independent of full-size examples.

The BSF files retain native column block ordering and a selected subgrid.
The SPC energy slice retains native rows; kx_ky uses a small orthogonal grid.
Full-resolution datasets remain in result_files/examples for numerical tests.
"""
from pathlib import Path

import numpy as np
import pytest

from .. import ResultFile

DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize(
    ("name", "keyword", "shape", "raw_size"),
    (
        ("SnTe_BLOCHSF_small.bsf", "BSF", (4, 4, 5), 240),
        ("SnTe_BLOCHSF_spol_small.bsf", "BSF-SPOL", (4, 4, 5), 320),
    ),
)
def test_compact_bsf_parses_all_native_blocks(name, keyword, shape, raw_size):
    result = ResultFile.from_file(DATA / name, unknown=False)
    assert result.KEYWORD() == keyword
    assert result.MODE() == "EK-REL"
    assert result.NE() == 4
    assert result.NK() == 5
    assert result.INDKDIR()[-1] == 5
    assert result.RAW_DATA().size == raw_size
    assert result.I().shape == shape
    channels = ("I_UP", "I_DOWN") if keyword == "BSF" else ("I_X", "I_Y", "I_Z")
    for channel in channels:
        assert result[channel]().shape == shape
        assert np.isfinite(result[channel]()).all()


@pytest.mark.parametrize(
    ("name", "mode", "np_points", "nt_points", "ne_points", "columns"),
    (
        ("Mo2Se4_ARPES_small.spc", "energy", 1, 1, 7, 8),
        ("Mo2Se4_ARPES_kx_ky_small.spc", "kx_ky", 5, 7, 1, 7),
    ),
)
def test_compact_arpes_parses_consistent_native_grid(
    name, mode, np_points, nt_points, ne_points, columns,
):
    result = ResultFile.from_file(DATA / name, unknown=False)
    assert result.KEYWORD() == "ARPES"
    assert result.MODE() == mode
    assert result.NP() == np_points
    assert result.NT() == nt_points
    assert result.NE() == ne_points
    assert result.RAW_DATA().shape == (np_points * nt_points * ne_points, columns)
    assert np.isfinite(result.TOTAL()).all()
    assert np.isfinite(result.POLARIZATION()).all()


@pytest.mark.parametrize(
    ("name", "mode"),
    (
        ("SnTe_BLOCHSF_small.bsf", "BSF"),
        ("SnTe_BLOCHSF_spol_small.bsf", "BSF-SPOL"),
    ),
)
def test_compactor_round_trip_bsf(tmp_path, name, mode):
    from .compact_spectral_data import compact_bsf

    reduced = compact_bsf((DATA / name).read_text(), max_e=3, max_k=3)
    path = tmp_path / name
    path.write_text(reduced)
    result = ResultFile.from_file(path, unknown=False)
    assert result.KEYWORD() == mode
    assert result.NE() == 3
    assert result.NK() == 3
    assert result.INDKDIR()[-1] == 3
    assert result.I().shape == (4, 3, 3)


@pytest.mark.parametrize(
    ("name", "mode"),
    (
        ("Mo2Se4_ARPES_small.spc", "energy"),
        ("Mo2Se4_ARPES_kx_ky_small.spc", "kx_ky"),
    ),
)
def test_compactor_round_trip_spc(tmp_path, name, mode):
    from .compact_spectral_data import compact_spc

    reduced = compact_spc((DATA / name).read_text(), max_e=4, max_k=3)
    path = tmp_path / name
    path.write_text(reduced)
    result = ResultFile.from_file(path, unknown=False)
    assert result.MODE() == mode
    assert np.isfinite(result.TOTAL()).all()
    if mode == "energy":
        assert result.NE() == 4
    else:
        assert result.NP() == 3
        assert result.NT() == 3


@pytest.mark.slow
@pytest.mark.parametrize(
    ("original", "expected_mode"),
    (
        ("FePt.new_BSFEK_BLOCHSF.bsf", "EK-REL"),
        ("FePt.new_BSFEK_BLOCHSF_spol.bsf", "EK-REL"),
        ("WSe2_BLOCHSF_spol_kx_ky.bsf", "CONST-E"),
        ("WSe2_ARPES.spc", "energy"),
        ("Au_kx_ky_ARPES_data.spc", "kx_ky"),
    ),
)
def test_compact_large_native_example(tmp_path, original, expected_mode):
    """Full-size source -> strict reducer -> standard reader, including channel geometry."""
    from .compact_spectral_data import compact_bsf, compact_spc

    source = Path(__file__).resolve().parents[1] / "examples" / original
    text = source.read_text()
    reduced = (compact_bsf(text, max_e=6, max_k=8)
               if source.suffix == ".bsf"
               else compact_spc(text, max_e=6, max_k=8))
    assert len(reduced) < len(text) / 4
    dest = tmp_path / original
    dest.write_text(reduced)
    result = ResultFile.from_file(dest, unknown=False)
    if source.suffix == ".bsf":
        assert result.MODE() == expected_mode
        if original.startswith("FePt.new_BSFEK_"):
            import re
            original_boundaries = re.findall(
                r"(?m)^\s*INDKDIR\s+(\d+)\s*$", text)
            compact_boundaries = re.findall(
                r"(?m)^\s*INDKDIR\s+(\d+)\s*$", reduced)
            assert len(compact_boundaries) == len(original_boundaries)
            assert list(map(int, compact_boundaries)) == sorted(
                set(map(int, compact_boundaries)))
            assert int(compact_boundaries[-1]) == result.NK()
        assert result.I().size > 0
        assert np.isfinite(result.I()).all()
        if result.KEYWORD() == "BSF":
            channels = ("I_UP", "I_DOWN")
        else:
            channels = ("I_X", "I_Y", "I_Z")
        for channel in channels:
            assert result[channel]().shape == result.I().shape
    else:
        assert result.MODE() == expected_mode
        assert result.TOTAL().size > 0
        assert np.isfinite(result.TOTAL()).all()
