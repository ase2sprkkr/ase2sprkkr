"""Structured result-file discovery must work without plotting extras."""
import subprocess
import sys


def test_discovery_does_not_import_matplotlib_or_scipy_interpolate():
    script = """
import builtins
real_import = builtins.__import__

def deny_optional(name, *args, **kwargs):
    if name == 'matplotlib' or name.startswith('matplotlib.'):
        raise ModuleNotFoundError("No module named 'matplotlib'", name='matplotlib')
    if name == 'scipy.interpolate' or name.startswith('scipy.interpolate.'):
        raise ModuleNotFoundError("No module named 'scipy.interpolate'", name='scipy.interpolate')
    return real_import(name, *args, **kwargs)

builtins.__import__ = deny_optional
from ase2sprkkr.result_files import ResultFile
definitions = ResultFile.definitions
assert 'spc' in definitions
assert 'bsf' in definitions
assert 'dos' in definitions
"""
    completed = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_only_arpes_spc_has_nexus_export_action():
    from ase2sprkkr.result_files.definitions.spc import ARPESResultFile
    from ase2sprkkr.result_files.definitions.band_spc import BandSPCResultFile
    from ase2sprkkr.result_files.definitions.spleed_spc import SPLEEDSPCResultFile

    assert "export_nexus" in ARPESResultFile.additional_actions
    assert "export_nexus" not in BandSPCResultFile.additional_actions
    assert "export_nexus" not in SPLEEDSPCResultFile.additional_actions
