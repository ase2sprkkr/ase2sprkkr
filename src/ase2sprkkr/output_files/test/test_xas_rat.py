"""Test that RATOutputFile.generate_data() reproduces the reference Fortran
XAS output for the Cu L2,3 example.
"""

import re
import numpy as np
import pytest
from scipy.interpolate import interp1d
import os

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from ase2sprkkr.output_files.test.init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DIRNAME = os.path.dirname(__file__)
RAT_FILE = os.path.join(DIRNAME, "..", "examples", "Cu_XAS_Cu_L23.rat")
REF_XAS_AGR = os.path.join(DIRNAME, "xas_Cu_L23_x.agr")
REF_SUMR_AGR = os.path.join(DIRNAME, "xas_Cu_L23_sumrule.agr")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def parse_agr(path):
    """Return dict {target_label: np.array shape (N,2)} for a Grace .agr file."""
    data = {}
    current = None
    buf = []
    with open(path) as f:
        for line in f:
            m = re.match(r"@target\s+(G\d+\.S\d+)", line)
            if m:
                if current and buf:
                    data[current] = np.array(buf)
                current = m.group(1)
                buf = []
                continue
            if current:
                if line.strip() == "&":
                    data[current] = np.array(buf)
                    current = None
                    buf = []
                elif not line.startswith("@") and not line.startswith("#"):
                    parts = line.split()
                    if len(parts) == 2:
                        buf.append([float(parts[0]), float(parts[1])])
    # Save last block if file did not terminate with '&'
    if current and buf:
        data[current] = np.array(buf)
    return data


@pytest.fixture(scope="module")
def generated_data():
    from ase2sprkkr.output_files.definitions.rat import RATOutputFile

    of = RATOutputFile.from_file(RAT_FILE)
    # Fortran default: WGAUSS = 0.0 (no Gaussian broadening)
    # Fortran uses Fuggle-Inglesfield core-hole widths — request by name
    return of.generate_data(gauss_width=0.0, core_hole_width="fuggle-inglesfield")


@pytest.fixture(scope="module")
def ref_xas():
    d = parse_agr(REF_XAS_AGR)
    pol = d.get("G1.S0")
    diff = d.get("G0.S0")
    if pol is None:
        raise RuntimeError(f"Reference XAS file {REF_XAS_AGR} missing G1.S0 block")
    if diff is None:
        diff = np.zeros_like(pol)
    return pol, diff  # (polarization_avg, difference)


@pytest.fixture(scope="module")
def ref_sumrule():
    d = parse_agr(REF_SUMR_AGR)
    spin = d.get("G1.S0")
    orb = d.get("G0.S0")
    if spin is None or orb is None:
        raise RuntimeError(
            f"Reference sumrule file {REF_SUMR_AGR} missing expected blocks"
        )
    return spin, orb  # (spin SSR, orbital OSR) – normalised to max=100


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestEnergyGrid:
    """The generated energy grid must match the reference exactly."""

    def test_n_points(self, generated_data, ref_xas):
        ref_pol, _ = ref_xas
        assert len(generated_data["ENERGY"]) == len(ref_pol), (
            "Number of energy points differs from reference"
        )

    def test_energy_range(self, generated_data, ref_xas):
        ref_pol, _ = ref_xas
        e = generated_data["ENERGY"]
        assert abs(e[0] - ref_pol[0, 0]) < 0.01, "Start energy differs from reference"
        assert abs(e[-1] - ref_pol[-1, 0]) < 0.01, "End energy differs from reference"

    def test_energy_values(self, generated_data, ref_xas):
        ref_pol, _ = ref_xas
        e = generated_data["ENERGY"]
        max_err = np.max(np.abs(e - ref_pol[:, 0]))
        assert max_err < 0.01, f"Max energy deviation {max_err:.4f} eV exceeds 0.01 eV"


class TestPolarizationSpectrum:
    """Polarization-averaged XAS spectrum (G1.S0 in xas_Cu_L23_x.agr).

    Tolerances:
    - Shape (median ratio) within 5 % of 1.0
    - Point-wise RMS relative error below 10 %
    - Peak position within 0.5 eV
    """

    def test_peak_position(self, generated_data, ref_xas):
        ref_pol, _ = ref_xas
        e = generated_data["ENERGY"]
        py = generated_data["POLARIZATION"][0]
        ref = ref_pol[:, 1]
        py_peak_e = e[np.argmax(py)]
        ref_peak_e = ref_pol[np.argmax(ref), 0]
        assert abs(py_peak_e - ref_peak_e) < 0.5, (
            f"Peak position: Python {py_peak_e:.3f} eV, ref {ref_peak_e:.3f} eV"
        )

    def test_peak_amplitude(self, generated_data, ref_xas):
        ref_pol, _ = ref_xas
        py = generated_data["POLARIZATION"][0]
        ref = ref_pol[:, 1]
        ratio = py.max() / ref.max()
        assert 0.80 < ratio < 1.20, (
            f"Peak amplitude ratio {ratio:.3f} outside [0.80, 1.20]"
        )

    def test_pointwise_median_ratio(self, generated_data, ref_xas):
        """Median point-wise ratio of Python/Reference should be close to 1."""
        ref_pol, _ = ref_xas
        e = generated_data["ENERGY"]
        py = generated_data["POLARIZATION"][0]
        ref = ref_pol[:, 1]
        # evaluate Python at reference energy points
        f = interp1d(e, py, bounds_error=False, fill_value=0)
        py_at_ref = f(ref_pol[:, 0])
        mask = ref > 0.05  # only compare where spectrum is significant
        ratio = py_at_ref[mask] / ref[mask]
        median = np.median(ratio)
        assert 0.90 < median < 1.10, (
            f"Median Python/Reference ratio {median:.3f} outside [0.90, 1.10]"
        )

    def test_pointwise_rms(self, generated_data, ref_xas):
        """RMS relative error over significant region should be below 15 %."""
        ref_pol, _ = ref_xas
        e = generated_data["ENERGY"]
        py = generated_data["POLARIZATION"][0]
        ref = ref_pol[:, 1]
        f = interp1d(e, py, bounds_error=False, fill_value=0)
        py_at_ref = f(ref_pol[:, 0])
        mask = ref > 0.05
        rms = np.sqrt(np.mean(((py_at_ref[mask] - ref[mask]) / ref[mask]) ** 2))
        assert rms < 0.15, f"RMS relative error {rms:.3f} exceeds 0.15"

    def test_broad_shape(self, generated_data, ref_xas):
        """Check key shape features: trough below peak, secondary hump, tail."""
        ref_pol, _ = ref_xas
        e = generated_data["ENERGY"]
        py = generated_data["POLARIZATION"][0]
        f = interp1d(e, py, bounds_error=True)

        # Main peak > 1.0 Mb
        assert py.max() > 1.0, f"Peak {py.max():.3f} Mb < 1.0 Mb"

        # Trough at E ≈ 13-17 eV should be below 0.75 Mb
        trough_indices = (e > 13) & (e < 17)
        trough_val = py[trough_indices].min()
        assert trough_val < 0.75, f"Trough {trough_val:.3f} Mb >= 0.75 Mb"

        # Hump at E ≈ 20-25 eV should be above trough
        hump_indices = (e > 19) & (e < 26)
        hump_val = py[hump_indices].max()
        assert hump_val > trough_val, "Secondary hump not above trough"


class TestDifferenceSpectrum:
    """For a non-magnetic Cu calculation the difference spectrum must be
    effectively zero (floating-point noise level)."""

    def test_difference_is_zero(self, generated_data, ref_xas):
        _, ref_diff_arr = ref_xas
        e = generated_data["ENERGY"]
        py_diff = generated_data["DIFFERENCE"][0]
        py_max = np.max(np.abs(py_diff))
        ref_diff = ref_diff_arr[:, 1]
        ref_max = np.max(np.abs(ref_diff))
        # Both should be negligible compared to the XAS peak
        assert ref_max < 1e-6, (
            f"Reference difference is non-trivial: max={ref_max:.2e} Mb"
        )
        assert py_max < 1e-6, f"Python difference is non-trivial: max={py_max:.2e} Mb"


class TestSumRuleSpectra:
    """SPIN and ORBIT spectra (sum rules).

    For a non-magnetic Cu calculation the XMCD is zero, so SPIN and ORBIT
    must also be effectively zero (both in Python and in the reference).
    The reference sumrule.agr values are just floating-point noise at
    the ~1e-7 Mb level.
    """

    def test_spin_and_orbit_near_zero(self, generated_data, ref_sumrule):
        ref_spin_arr, ref_orb_arr = ref_sumrule
        if "SPIN" not in generated_data:
            pytest.skip("No SPIN data (n_ktypes <= 1)")

        e = generated_data["ENERGY"]
        py_spin = generated_data["SPIN"]
        py_orb = generated_data["ORBIT"]

        # Reference normalised to max=100; check that 100 corresponds to a physically negligible value
        ref_spin = ref_spin_arr[:, 1]
        ref_orb = ref_orb_arr[:, 1]
        py_scale = max(np.max(np.abs(py_spin)), np.max(np.abs(py_orb)))
        ref_scale = max(np.max(np.abs(ref_spin)), np.max(np.abs(ref_orb)))

        # Both should be negligible compared to the polarization-averaged XAS peak
        xas_peak = generated_data["POLARIZATION"][0].max()
        assert py_scale / xas_peak < 1e-5, (
            f"Python SPIN/ORBIT not near zero: {py_scale:.2e} vs XAS peak {xas_peak:.3f}"
        )
        assert ref_scale <= 100.1, (
            f"Reference sumrule unexpectedly large: {ref_scale:.1f}"
        )
