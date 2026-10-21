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
from ...gui.plot import change_default_kwargs, colormesh, lineplot, Multiplot


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
        if phi_range is not None and len(phi_range) != 2:
            raise ValueError("phi_range must contain (start_angle, end_angle)")
        phi = None if phi_range is None else np.linspace(
            phi_range[0], phi_range[1], self.NP()
        )

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

        live_options = None
        if multi and live:
            live_options = {}
            if self.NT() > 1:
                live_options["direction"] = {"options": ("NP", "NT")}
            live_options["slice"] = {
                "slider": lambda state: (
                    range(self.NT()) if state.get("direction") == "NT"
                    else indices
                ),
                "label": "Slice",
                "format": lambda value, state: (
                    f"{phi[value]:g}°"
                    if phi is not None and state.get("direction", "NP") == "NP"
                    else str(value)
                ),
            }

        with Multiplot(
            layout=layout,
            figsize=figsize,
            latex=latex,
            filename=filename,
            show=show,
            dpi=dpi,
            separate_plots=separate_plots,
            layout_kind=layout_kind,
            live=live,
            live_options=live_options,
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
        # Polarization is derived from the modified spin intensities.
        up, down = self.UP(), self.DOWN()
        total = up + down
        polarization = self.POLARIZATION()
        np.divide(up - down, total, out=polarization, where=total != 0)
        polarization[total == 0] = 0

    def _check_arithmetic(self, other):
        """Do not combine spectra with different sampled coordinates."""
        super()._check_arithmetic(other)
        if isinstance(other, numbers.Number):
            return
        if self.MODE() != other.MODE():
            raise ValueError("ARPES modes are different")
        if self.RAW_DATA().shape != other.RAW_DATA().shape:
            raise ValueError("ARPES data shapes are different")
        if self.NT() != other.NT() or self.NP() != other.NP():
            raise ValueError("ARPES grid dimensions are different")
        axes = ("ENERGY", "THETA", "K") if self.MODE() == "energy" else ("KX", "KY")
        for name in axes:
            left, right = self[name](), other[name]()
            if left.shape != right.shape or not np.allclose(left, right):
                raise ValueError(f"ARPES {name} coordinates are different")


class ARPESDefinition(ResultFileDefinition):
    result_class = ARPESResultFile


def create_definition():

    def plot(option, azimuth_index=None, azimuth_indices=None,
             live=False, phi_values=None, controls=None, **kwargs):
        c = option._container
        kw = {
            "show_zero_line": False,
            "mode": "from_zero",
            "norm": "log",
            "vmax": np.max([c.TOTAL().max(), c.UP().max(), c.DOWN().max()]),
        }
        if option.name == "POLARIZATION":
            # Polarization is signed and dimensionless, unlike intensities.
            kw.update(mode="centered", norm="lin", vmax=1.0)

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
        if c.MODE() == "energy" and c.NP() > 1:
            selected = tuple(azimuth_indices) if azimuth_indices is not None else (
                0 if azimuth_index is None else azimuth_index,
            )
            if controls is not None:
                theta, energy = x, y
                phi = (
                    np.arange(c.NP()) if phi_values is None
                    else np.asarray(phi_values)
                )

                def section(data, state):
                    if state.get("direction", "NP") == "NT":
                        return data[:, :, state["slice"]]
                    return data[:, state["slice"], :]

                x = lambda state: (
                    np.broadcast_to(phi, (energy.shape[0], c.NP()))
                    if state.get("direction", "NP") == "NT"
                    else section(theta, state)
                )
                y = lambda state: section(energy, state)
                raw_intensity = intensity
                intensity = lambda state: section(raw_intensity, state)

                def xlabel(state):
                    if state.get("direction", "NP") == "NT":
                        return r"$\phi$ (deg)" if phi_values is not None else "Azimuth index (NP)"
                    # A single polar-angle point is an energy spectrum;
                    # its horizontal axis must remain energy while sliding.
                    return kw["ylabel"] if c.NT() == 1 else r"$\theta$ (deg)"

                kw["xlabel"] = xlabel(controls.state)
                # Labels belong to the plot; changing direction changes x.
                axis = kwargs.get("axis")
                # A singleton mesh is rendered as a line: its horizontal
                # coordinate can be energy, not the 2D mesh x coordinate.
            else:
                x, y, intensity = (
                    data[:, selected[0], :] for data in (x, y, intensity)
                )
        kw.update(kwargs)
        first = intensity(controls.state) if callable(intensity) else intensity
        if 1 in first.shape:
            # A singleton mesh contains no Gouraud cells. Show the actual
            # one-dimensional spectrum rather than an empty pseudocolor plot.
            if c.MODE() == "energy":
                variable_is_y = first.shape[1] == 1
            else:
                initial_x = x(controls.state) if callable(x) else x
                initial_y = y(controls.state) if callable(y) else y
                variable_is_y = np.ptp(initial_y) > np.ptp(initial_x)
            coordinate = y if variable_is_y else x
            if callable(coordinate):
                data = coordinate
                coordinate = lambda index: data(index).ravel()
            else:
                coordinate = coordinate.ravel()
            if callable(intensity):
                data = intensity
                intensity = lambda index: data(index).ravel()
            else:
                intensity = intensity.ravel()
            kw["xlabel"] = kw["ylabel"] if variable_is_y else kw["xlabel"]
            kw["ylabel"] = option.name
            for key in ("show_zero_line", "mode", "norm", "vmax", "vmin", "colormap"):
                kw.pop(key, None)
            lineplot(coordinate, intensity, controls=controls, **kw)
        else:
            if controls is not None and callable(x) and kwargs.get("axis") is not None:
                axis = kwargs["axis"]
                controls.connect(lambda state: axis.set_xlabel(xlabel(state)))
            colormesh(x, y, intensity, controls=controls, **kw)

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
