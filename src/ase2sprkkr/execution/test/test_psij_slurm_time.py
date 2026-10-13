import re

from ase2sprkkr.execution.psij import _slurm_time


def test_slurm_accounting_time_drops_timezone_offset():
    value = _slurm_time("2026-10-13T08:00:00+02:00")

    assert re.fullmatch(r"2026-10-13T\d{2}:\d{2}:\d{2}", value)
    assert "+" not in value
