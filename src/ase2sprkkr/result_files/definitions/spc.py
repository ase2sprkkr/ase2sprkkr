from ..result_files_definitions import (
    ResultFileValueDefinition as V,
    create_result_file_definition,
    ResultFileDefinition,
)
from typing import Optional
import numpy as np

from ..result_files import Arithmetic, CommonResultFile
from ...common.grammar_types import NumpyArray, Prefixed
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV, GeneratedValueDefinition as GV
from ...common.configuration_definitions import switch
from ...gui.plot import change_default_kwargs, colormesh, Multiplot


class ARPESResultFile(CommonResultFile, Arithmetic):
    def to_nexus(self, filename, *, input_parameters=None, input_file=None,
                 metadata=None):
        """Export calculated energy/angle ARPES data and its provenance to NeXus."""
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
        **kwargs,
    ):
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
            mp.plot(self.TOTAL)
            mp.plot(self.UP)
            mp.plot(self.DOWN)
            mp.plot(self.POLARIZATION)

    def _arithmetic_values(self):
        return [("RAW_DATA", (slice(None), slice(2, 6 if self.MODE() == "energy" else 5)))]

    def _assert_arithmetic(self, other):
        """Check, that the file can be summed/subtracked from an other file"""
        assert self.MODE() == other.MODE()
        if self.MODE() == "energy":
            assert np.allclose(other.ENERGY(), self.ENERGY())
            assert np.allclose(other.THETA(), self.THETA())
        else:
            assert np.allclose(other.KX(), self.KX())
            assert np.allclose(other.KY(), self.KY())


class ARPESDefinition(ResultFileDefinition):
    result_class = ARPESResultFile


def create_definition():

    def plot(option, **kwargs):
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

        kw.update(kwargs)
        colormesh(x, y, option(), **kw)

    def i(j):
        return slice(None), j

    def determinant(data: np.ndarray, section: object) -> np.ndarray:
        """Reshape determinant ``data`` using dimensions from ``section``."""
        if section["MODE"]() == "energy":
            return data[:, 7].reshape(section["NE"](), section["NT"]())
        return data[:, 6].reshape(section["NP"](), section["NT"]())

    def ii(j):
        def fn(data: np.ndarray, section: object) -> np.ndarray:
            """Select column ``j`` from ``data`` and shape it from ``section``."""
            first = section["NE"]() if section["MODE"]() == "energy" else section["NP"]()
            return data[:, j].reshape(first, section["NT"]())

        return fn

    def mode(section: object, _key: object = None) -> str:
        """Derive the output mode from ``section``."""
        return "energy" if section["NE"]() > 1 else "kx_ky"

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
                        NV("THETA", "RAW_DATA", i(0), ("NE", "NT"), info="Photoelectron polar angle in degrees"),
                        NV("ENERGY", "RAW_DATA", i(1), ("NE", "NT"), info="Initial-state energy E - E_F (eV)"),
                        NV(
                            "K", "RAW_DATA", i(6), ("NE", "NT"),
                            info="Native parallel momentum coordinate (1/angstrom; signed angular scan or K-vector displacement)",
                        ),
                    ],
                    "kx_ky": [
                        NV("KX", "RAW_DATA", i(0), ("NP", "NT"), info="k_x coordinate (1/Å)"),
                        NV("KY", "RAW_DATA", i(1), ("NP", "NT"), info="k_y coordinate (1/Å)"),
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
