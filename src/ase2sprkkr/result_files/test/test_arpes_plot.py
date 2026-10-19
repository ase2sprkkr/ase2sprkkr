"""Quick real-SPC regression tests for ARPES figure content and sliders."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from ase2sprkkr.result_files import ResultFile


_EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "Mo2Se4_ARPES_ARPES_data.spc"


def _artists(figure):
    controls = getattr(figure, "_ase2sprkkr_controls", None)
    control_axes = set() if controls is None else {
        widget.ax for widget in controls.widgets.values()
    }
    return [
        axis for axis in figure.axes
        if axis not in control_axes and (axis.lines or axis.collections)
    ]


def test_real_single_angle_spectrum_is_not_blank():
    spc = ResultFile.from_file(str(_EXAMPLE), unknown=False)
    assert spc.NP() == spc.NT() == 1
    try:
        spc.plot(show=False)
        figure = plt.gcf()
        panels = _artists(figure)
        assert len(panels) == 4
        assert all(len(ax.lines) == 1 for ax in panels)
        for axis in panels:
            line = axis.lines[0]
            assert len(line.get_xdata()) == spc.NE()
            assert len(line.get_ydata()) == spc.NE()
            assert np.max(line.get_ydata()) > 0
        assert len(figure.axes) == 4  # NP=1: no artificial slider
    finally:
        plt.close("all")


def _multi_azimuth_file(tmp_path):
    """Turn a small real NP=1 SPC fixture into an NP=3 native SPC input."""
    text = _EXAMPLE.read_text()
    header, payload = text.split("#######################", 1)
    header = header.replace("NP                 1", "NP                 3")
    assert "NP                 3" in header
    output = tmp_path / "synthetic_np3.spc"
    lines = []
    for line in payload.splitlines():
        if not line.strip():
            continue
        values = np.fromstring(line, sep=" ")
        assert len(values) == 8
        for ip in range(3):
            modified = values.copy()
            modified[2:5] *= ip + 1
            modified[5] = (
                (modified[3] - modified[4]) /
                (modified[3] + modified[4])
            )
            lines.append(" ".join(f"{x:.9E}" for x in modified))
    output.write_text(header + "#######################\n" + "\n".join(lines) + "\n")
    return output


def test_real_parser_multi_azimuth_slider_draws_all_four_panels(tmp_path):
    spc = ResultFile.from_file(str(_multi_azimuth_file(tmp_path)), unknown=False)
    assert spc.NP() == 3
    assert spc.TOTAL().shape == (spc.NE(), 3, 1)
    try:
        spc.plot(show=False, live=True)
        fig = plt.gcf()
        assert len(_artists(fig)) == 4
        controls = fig._ase2sprkkr_controls
        assert controls._slider_values["slice"] == (0, 1, 2)
        lines = [ax.lines[0] for ax in _artists(fig)]
        before = np.asarray(lines[0].get_ydata()).copy()
        controls.widgets["slice"].set_val(2)
        after = np.asarray(lines[0].get_ydata())
        np.testing.assert_allclose(after, 3 * before)
    finally:
        plt.close("all")


def test_multi_azimuth_static_slice_count(tmp_path):
    spc = ResultFile.from_file(str(_multi_azimuth_file(tmp_path)), unknown=False)
    try:
        spc.plot(show=False, live=False, azimuths=2)
        fig = plt.gcf()
        assert len(_artists(fig)) == 8
        assert not hasattr(fig, "_ase2sprkkr_controls")
    finally:
        plt.close("all")



def test_controls_switch_np_nt_on_real_parser_grid(tmp_path):
    """A parsed NP x NT spectrum changes slice direction and slider limits."""
    header, payload = _EXAMPLE.read_text().split("#######################", 1)
    header = header.replace("NT                 1", "NT                 3")
    header = header.replace("NP                 1", "NP                 4")
    lines = []
    for line in payload.splitlines():
        if not line.strip():
            continue
        values = np.fromstring(line, sep=" ")
        for ip in range(4):
            for it in range(3):
                point = values.copy()
                point[0] = -10 + 10 * it
                point[2:5] *= (ip + 1) * (it + 1)
                lines.append(" ".join(f"{v:.9E}" for v in point))
    path = tmp_path / "two_angles.spc"
    path.write_text(header + "#######################\n" + "\n".join(lines) + "\n")
    spc = ResultFile.from_file(str(path), unknown=False)
    assert spc.TOTAL().shape == (spc.NE(), 4, 3)
    try:
        spc.plot(show=False, live=True)
        fig = plt.gcf()
        controls = fig._ase2sprkkr_controls
        assert len(_artists(fig)) == 4
        assert controls.state == {"direction": "NP", "slice": 0}
        assert len(controls._slider_values["slice"]) == 4
        assert _artists(fig)[0].get_xlabel() == r"$\theta$ (deg)"
        controls.widgets["slice"].set_val(3)
        controls.widgets["direction"].set_active(1)
        assert controls.state == {"direction": "NT", "slice": 0}
        assert len(controls._slider_values["slice"]) == 3
        assert _artists(fig)[0].get_xlabel() == "Azimuth index (NP)"
        controls.widgets["slice"].set_val(2)
        assert controls.state["slice"] == 2
        assert len(_artists(fig)) == 4
    finally:
        plt.close("all")
