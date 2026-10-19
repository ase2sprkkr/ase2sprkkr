from ..result_files_definitions import (
    ResultFileValueDefinition as V,
    create_result_file_definition,
    ResultFileDefinition,
)
from typing import Optional
import numbers
import numpy as np

from ..result_files import Arithmetic, CommonResultFile
from ...common.grammar_types import NumpyArray, Prefixed
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV, GeneratedValueDefinition as GV
from ...common.configuration_definitions import switch
from ...gui.plot import change_default_kwargs, colormesh, Multiplot


# Native SPC columns are shared by the reader and the NeXus writer.
SPC_COLUMNS = {
    "energy": ("THETA", "ENERGY", "TOTAL", "UP", "DOWN",
               "POLARIZATION", "K", "DETERMINANT"),
    "kx_ky": ("KX", "KY", "TOTAL", "UP", "DOWN",
              "POLARIZATION", "DETERMINANT"),
}


class ARPESResultFile(CommonResultFile, Arithmetic):
    def to_nexus(self, filename, *, input_parameters=None, input_file=None,
                 metadata=None):
        """Export an ARPES spectrum together with its available input data."""
        from ..nexus import write_arpes_nexus

        return write_arpes_nexus(
            self, filename, input_parameters=input_parameters,
            input_file=input_file, metadata=metadata
        )

    def plot(
        self,
        layout=(2, 2),
        figsize=(10, 6),
        latex=None,
        filename: Optional[str] = None,
        show: Optional[bool] = None,
        dpi=800,
        separate_plots=False,
        layout_kind="constrained",
        azimuths=None,
        phi_range=None,
        input_parameters=None,
        live=None,
        **kwargs,
    ):
        """Plot energy-angle slices, optionally choosing azimuths.

        azimuths: integer number of evenly distributed slices, or a sequence
        of phi angles in degrees. Physical angles need phi_range or the
        corresponding input parameters. Without a selection, static plots
        use the first slice; live plots can switch between all slices.
        """
        multi = self.MODE() == "energy" and self.NP() > 1
        if live is None:
            live = filename is None and show is not False
        if live and separate_plots:
            raise ValueError("A shared live slider requires one figure")

        if input_parameters is not None and phi_range is None:
            phi_range = input_parameters.SPEC_EL.PHI()
        phi = None if phi_range is None else np.linspace(*phi_range, self.NP())

        if not multi:
            if azimuths is not None:
                raise ValueError("azimuths selection requires NP > 1")
            indices = (None,)
            live = False
        elif azimuths is None:
            indices = tuple(range(self.NP())) if live else (0,)
        elif isinstance(azimuths, numbers.Integral) and not isinstance(azimuths, bool):
            if azimuths <= 0:
                raise ValueError("azimuths must be a positive number of slices")
            indices = tuple(np.rint(np.linspace(
                0, self.NP() - 1, min(azimuths, self.NP())
            )).astype(int))
        else:
            if phi is None:
                raise ValueError("Azimuth angles require phi_range or input_parameters")
            indices = tuple(int(np.abs(phi - angle).argmin()) for angle in azimuths)
            if not indices:
                raise ValueError("azimuths must not be empty")

        # Both paths use the same plotting operation. Live mode switches
        # among slices; static mode renders the requested slices as panels.
        slices = (indices,) if live else tuple((index,) for index in indices)
        if multi and not live and layout == (2, 2) and len(slices) > 1:
            layout = (2 * len(slices), 2)
        if layout[0] * layout[1] < 4 * len(slices):
            raise ValueError("The layout has too few axes for the selected slices")

        with Multiplot(
            layout=layout,
            figsize=figsize,
            latex=latex,
            filename=filename,
            show=show,
            dpi=dpi,
            separate_plots=separate_plots,
            layout_kind=layout_kind,
            **kwargs,
        ) as mp:
            for selected in slices:
                for option in (self.TOTAL, self.UP, self.DOWN, self.POLARIZATION):
                    with mp.new_axis(name=option.name) as (axis, plot_kwargs):
                        option.plot(
                            azimuth_indices=selected if multi else None,
                            live=live, phi_values=phi, axis=axis, **plot_kwargs
                        )

    def spc_columns(self):
        """Native column ordering for this file's detected SPC mode."""
        return SPC_COLUMNS[self.MODE()]

    def _arithmetic_values(self):
        # Spin polarization is a ratio, not another additive intensity.
        return [("RAW_DATA", (slice(None), slice(2, 5)))]

    def _do_arithmetic(self, func, other):
        super()._do_arithmetic(func, other)
        if func in ("__iadd__", "__isub__"):
            up, down = self.UP(), self.DOWN()
            np.divide(
                up - down, up + down,
                out=self.POLARIZATION(), where=(up + down) != 0
            )
            self.POLARIZATION()[(up + down) == 0] = 0

    def _assert_arithmetic(self, other):
        """Check, that the file can be summed/subtracked from an other file"""
        assert self.MODE() == other.MODE()
        assert self.RAW_DATA().shape == other.RAW_DATA().shape
        if self.MODE() == "energy":
            assert np.allclose(other.ENERGY(), self.ENERGY())
            assert np.allclose(other.THETA(), self.THETA())
        else:
            assert np.allclose(other.KX(), self.KX())
            assert np.allclose(other.KY(), self.KY())


class ARPESDefinition(ResultFileDefinition):
    result_class = ARPESResultFile


def create_definition():

    def plot(option, azimuth_index=None, azimuth_indices=None,
             live=False, phi_values=None, **kwargs):
        c = option._container
        kw = {
            "show_zero_line": False,
            "mode": "from_zero",
            "norm": "log",
            "vmax": np.max([c.TOTAL().max(), c.UP().max(), c.DOWN().max()]),
        }

        if c.MODE() == "energy":
            x = c.THETA()
            y = c.ENERGY()
            kw["xlabel"] = r"$\theta$ (deg)"
            kw["ylabel"] = r"$E-E_{\rm F}$ (eV)"
        else:
            x = c.KX()
            y = c.KY()
            kw["xlabel"] = r"$k_x$ (Å$^{-1}$)"
            kw["ylabel"] = r"$k_y$ (Å$^{-1}$)"

        intensity = option()
        selected = None
        labels = None
        if c.MODE() == "energy" and c.NP() > 1:
            selected = tuple(azimuth_indices) if azimuth_indices is not None else (
                0 if azimuth_index is None else azimuth_index,
            )
            x, y, intensity = (
                lambda index, data=data: data[:, index, :]
                for data in (x, y, intensity)
            )
            if phi_values is not None:
                labels = tuple(f"{phi_values[i]:g}°" for i in selected)
        kw.update(kwargs)
        colormesh(
            x, y, intensity, live=live, slider_values=selected,
            slider_labels=labels, slider_label="Azimuth", **kw
        )

    def ii(j):
        def fn(data: np.ndarray, section: object) -> np.ndarray:
            if section["MODE"]() == "energy":
                shape = (section["NE"](), section["NT"]())
                if section["NP"]() > 1:
                    shape = (section["NE"](), section["NP"](), section["NT"]())
            else:
                shape = (section["NP"](), section["NT"]())
            return data[:, j].reshape(shape)
        return fn

    def determinant(data: np.ndarray, section: object) -> np.ndarray:
        return ii(7 if section["MODE"]() == "energy" else 6)(data, section)

    def mode(section: object, _key: object = None) -> str:
        """Detect native mode by record width (not by the number of energies)."""
        columns = section["RAW_DATA"]().shape[1]
        for mode_name, names in SPC_COLUMNS.items():
            if len(names) == columns:
                return mode_name
        raise ValueError(f"Unrecognized ARPES SPC record width: {columns}")

    definition = create_result_file_definition(
        "ARPES",
        [
            V("NT", int),
            V("NP", int),
            V("COMMENT", Prefixed("#"), name_in_grammar=False),
            V("RAW_DATA", NumpyArray(), name_in_grammar=False),
            GV("MODE", mode),
            *switch(
                "MODE",
                {
                    "energy": [
                        NV("THETA", "RAW_DATA", ii(0), info="Photoelectron polar angle in degrees"),
                        NV("ENERGY", "RAW_DATA", ii(1), info="Initial-state energy E - E_F (eV)"),
                        NV(
                            "K", "RAW_DATA", ii(6),
                            info="Native parallel momentum coordinate (1/angstrom)",
                        ),
                    ],
                    "kx_ky": [
                        NV("KX", "RAW_DATA", ii(0), info="k_x coordinate (1/Å)"),
                        NV("KY", "RAW_DATA", ii(1), info="k_y coordinate (1/Å)"),
                    ],
                },
            ),
            NV(
                "TOTAL",
                "RAW_DATA",
                ii(2),
                info="Total intensity",
                plot=change_default_kwargs(plot, title=r"Total intensity", colormap="gray"),
            ),
            NV(
                "UP",
                "RAW_DATA",
                ii(3),
                info="Spin up",
                plot=change_default_kwargs(plot, title=r"Spin up", colormap="gray"),
            ),
            NV(
                "DOWN",
                "RAW_DATA",
                ii(4),
                info="Spin down",
                plot=change_default_kwargs(plot, title=r"Spin down", colormap="gray"),
            ),
            NV(
                "POLARIZATION",
                "RAW_DATA",
                ii(5),
                info="Spin polarization",
                plot=change_default_kwargs(plot, colormap="bwr", title=r"Spin polarization"),
            ),
            NV("DETERMINANT", "RAW_DATA", determinant),
        ],
        cls=ARPESDefinition,
        info="ARPES (Angle-resolved photoemission spectroscopy) output (.spc file).",
    )
    return definition


definition = create_definition()
