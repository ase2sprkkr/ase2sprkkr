from ..output_files_definitions import (
    OutputFileValueDefinition as V,
    create_output_file_definition,
    OutputFileDefinition,
    Separator,
)
from typing import Optional
import numpy as np

from ..output_files import Arithmetic, CommonOutputFile
from ...common.grammar_types import Array, NumpyArray, Keyword
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV
from ...common.configuration_definitions import gather, switch
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from functools import lru_cache
from ase.units import Rydberg
from ...gui.plot import colormesh, Multiplot


class BSFOutputFile(CommonOutputFile, Arithmetic):
    plot_parameters = {"sites", "fermi", "seperate_plots"}

    """
    Output file for Bloch spectral functions
    """

    def plot(
        self,
        layout=(2, 2),
        figsize=(10, 6),
        latex=None,
        filename: Optional[str] = None,
        show: Optional[bool] = None,
        dpi=600,
        layer=None,
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
            if self.KEYWORD() in ["BSF-SPN", "BSF-SPOL"]:
                mp.plot(self.I)
                mp.plot(self.I_X)
                mp.plot(self.I_Y)
                mp.plot(self.I_Z)
            else:
                mp.plot(self.I)
                mp.plot(self.I_UP)
                mp.plot(self.I_DOWN)

    _arithmetic_values = [("RAW_DATA", slice(None))]

    def _assert_arithmetic(self, other):
        """Check, that the file can be summed/subtracked from an other file"""
        assert self.MODE() == other.MODE()
        assert self.KEYWORD() == other.KEYWORD()
        for i in ["SHAPE", "NQ_EFF"]:
            assert self[i]() == other[i]()
        if self.MODE() == "EK_REL":
            assert self.NK() == other.NK() and self.NE() == other.NE()
            assert np.allclose(self.K(), other.K())
        else:
            assert self.NK1() == other.NK1() and self.NK2() == other.NK2()
            assert np.allclose(self.VECK1(), other.VECK1())
            assert np.allclose(self.VECK2(), other.VECK2())


class BSFDefinition(OutputFileDefinition):
    result_class = BSFOutputFile


def create_definition():
    cmap1 = plt.cm.bone_r(np.linspace(0.0, 1.0, 256))
    cmap2 = plt.cm.hot(np.linspace(0.0, 1.0, 256))
    cmap = np.vstack((cmap1, cmap2))
    mymap = mcolors.LinearSegmentedColormap.from_list("my_colormap", cmap)

    def plot(title, colormap="bwr", negative=True):

        def plot(option, colormap=colormap, sites=None, fermi=None, **kwargs):
            c = option._container
            mesh = c.MESH()
            data = option()
            if fermi is True:
                fermi = 0.5

            def check_sites(site):
                limit = data.shape[0]
                if site < 0 or site >= limit:
                    raise ValueError(
                        f"Site number should be between {0} and {limit - 1}."
                    )

            if isinstance(sites, int):
                check_sites(sites)
                data = data[sites]
            else:
                if isinstance(sites, list):
                    for i in sites:
                        check_sites(i)
                elif isinstance(sites, slice):
                    check_sites(sites.start)
                    check_sites(sites.stop)
                    if sites.start >= sites.stop:
                        raise ValueError("Empty sites range.")
                elif not sites is None:
                    raise ValueError(f"I can not filter sites by {sites}.")

                data = data[sites].sum(axis=0)

            if negative:
                vmax = max(np.abs(np.max(data)), np.abs(np.min(data)))
                vmin = -vmax
            else:
                vmin = 0
                vmax = np.max(data)
            k = mesh[0, 0]
            kw = {
                "vmin": vmin,
                "vmax": vmax,
                "colormap": colormap,
                "colorbar": True,
                "title": title,
                "show_zero_line": fermi,
            }

            if c.MODE() == "CONST-E":

                def ticks(v):
                    return np.array2string(np.array(v), precision=3, separator=",")

                kw.update(
                    {
                        "xticks": [0.0, 1.0],
                        "xticklabels": [
                            "" if c.VECK_START() is None else v(c.VECK_START()),
                            ticks(c.VECK1()),
                        ],
                        "xlabel": r"Kx",
                        "ylabel": r"Ky",
                        "yticks": [1.0],
                        "yticklabels": [ticks(c.VECK2())],
                    }
                )

                def callback(ax):
                    ax.tick_params(axis="y", rotation=90)
                    ax.tick_params(axis="both", labelsize=8)
            else:
                kw.update(
                    {
                        "ylabel": r"$E-E_{\rm F}$ (eV)",
                        "xticks": np.insert(k[[x - 1 for x in c.INDKDIR()]], 0, 0),
                        "xticklabels": [],
                        "xlabel": r"K",
                    }
                )

                def callback(ax):
                    for index in c.INDKDIR()[:-1]:
                        ax.plot(
                            [k[index - 1], k[index - 1]],
                            [mesh[1, 0, 1], mesh[1, -1, 1]],
                            color="black",
                            lw=0.5,
                        )

            kw.update(kwargs)
            colormesh(*mesh, data, callback=callback, **kw)

        return plot

    reorder = (1, 0, 2)

    @lru_cache(maxsize=2)
    def k_points(start, end, num):
        return np.linspace(start, end, num)

    @lru_cache(maxsize=2)
    def energy_points(start, end, fermi, num):
        return np.linspace(start - fermi, end - fermi, num) * Rydberg

    def i(type):
        def index(data, c):
            """
            Returns data in the shape
            ('NE/NK1','NQ_EFF', 'NK(2)')
            for a given type.
            Reorder parameter then change the order of the axes to
            ('NQ_EFF', 'NE/NK1', 'NK(2)')

            Data structure:
              IX,Y,Z (BSF-SPOL/SPN)
                  k-e: NE, type, NQ, NK  types(I,x,y,z)
                  k-k: type, K1, NQ, K2, types(I,x,y,z)

             Iup,dn (BSF) two sets
                  k-e: NE, type, NQ, NK   types(u,d)
                       NE, type, NQ, NK   types(I)
                  k-k: NK, type, NQ, NK   types(u,d)
                       NK, NQ, NK   types(I)
            """
            nq = c.NQ_EFF()

            ekrel = c.MODE() == "EK-REL"
            nk2 = c.NK() if ekrel else c.NK2()
            if c.KEYWORD() == "BSF":
                nk1 = c.NE() if ekrel else c.NK1()
                limit = nq * nk1 * nk2 * 2
                if type >= 0:
                    return data[:limit].reshape(nk1, 2, nq, nk2)[:, type]
                else:
                    return data[limit:].reshape(nk1, nq, nk2)
            else:
                if ekrel:
                    return data.reshape(c.NE(), 4, nq, nk2)[:, type + 1]
                else:
                    return data.reshape(4, c.NK1(), nq, nk2)[type + 1]

        return index

    definition = create_output_file_definition(
        Keyword("BSF-SPOL", "BSF-SPN", "BSF"),
        [
            Separator(),
            V("DATASET", str, written_name="#DATASET"),
            V("MODE", Keyword("EK-REL", "CONST-E")),
            *switch(
                "MODE",
                {
                    "EK-REL": [
                        V(
                            "NE_a",
                            int,
                            written_name="NE",
                            info="Number of energies (the second axis)",
                            is_hidden=True,
                        ),
                        V("NK", int, info="Number of K points (the same as NK2)"),
                        Separator(),
                        *gather(V("EMIN", float), V("EMAX", float)),
                        Separator(),
                        V("NKDIR", int),
                        V(
                            "LBLKDIR",
                            NumpyArray(
                                written_shape=(-1, 1),
                                shape=(-1,),
                                lines="NKDIR",
                                dtype="line",
                            ),
                            name_in_grammar=False,
                        ),
                        Separator(),
                        V("INDKDIR", int, is_repeated=True),
                        Separator(),
                        V(
                            "NK_a",
                            int,
                            written_name="NK",
                            info="Number of K points (the last axis)",
                            is_hidden=True,
                        ),
                        V(
                            "K",
                            NumpyArray(written_shape=(-1, 1), shape=(-1,), lines="NK"),
                            name_in_grammar=False,
                        ),
                        Separator(),
                        V(
                            "E",
                            Array(float),
                            default_value_from_container=lambda o: energy_points(
                                o.EMIN(), o.EMAX(), o.EFERMI(), o.NE()
                            ),
                            is_stored=False,
                            info="Energy (relative to Fermi energy)",
                        ),
                    ],
                    "CONST-E": [
                        V("NK1", int, info="Number of K points (the first axis)"),
                        V("NK2", int, info="Number of K points (the second axis)"),
                        V("ERYD", Array(float, length=2)),
                        Separator(),
                        V(
                            "VECK_START",
                            Array(float, length=3),
                            is_stored=False,
                            is_optional=True,
                        ),
                        V("NK1_a", int, written_name="NK1", is_hidden=True),
                        V("VECK1", Array(float, length=3)),
                        V("NK2_a", int, written_name="NK2", is_hidden=True),
                        V("VECK2", Array(float, length=3)),
                        Separator(),
                        V(
                            "K1",
                            Array(float),
                            init_by_default=True,
                            default_value_from_container=lambda o: k_points(
                                0.0,
                                1.0,
                                o.NK1(),
                            ),
                            is_stored=False,
                            info="First axis for the data",
                        ),
                        V(
                            "K2",
                            Array(float),
                            init_by_default=True,
                            default_value_from_container=lambda o: k_points(
                                0.0,
                                1.0,
                                o.NK2(),
                            ),
                            is_stored=False,
                            info="Second axis for the data",
                        ),
                    ],
                },
            ),
            V(
                "MESH",
                NumpyArray(float),
                is_stored=False,
                init_by_default=True,
                default_value_from_container=lambda c: (
                    np.meshgrid(c.K1(), c.K2())
                    if c.MODE() == "CONST-E"
                    else np.meshgrid(c.K(), c.E())
                ),
            ),
            V(
                "RAW_DATA",
                NumpyArray(written_shape=(-1, 1), shape=(-1,)),
                name_in_grammar=False,
            ),
            *switch(
                "KEYWORD",
                {
                    "BSF": [
                        NV(
                            "I_UP",
                            "RAW_DATA",
                            i(0),
                            reorder=reorder,
                            plot=plot(title="Spin up", negative=True, colormap="Reds"),
                        ),
                        NV(
                            "I_DOWN",
                            "RAW_DATA",
                            i(1),
                            reorder=reorder,
                            plot=plot(
                                title="Spin down", negative=True, colormap="Blues"
                            ),
                        ),
                    ],
                    "BSF-SPOL": [
                        NV(
                            "I_X",
                            "RAW_DATA",
                            i(0),
                            reorder=reorder,
                            plot=plot(title=r"$\sigma_x$"),
                        ),
                        NV(
                            "I_Y",
                            "RAW_DATA",
                            i(1),
                            reorder=reorder,
                            plot=plot(title=r"$\sigma_y$"),
                        ),
                        NV(
                            "I_Z",
                            "RAW_DATA",
                            i(2),
                            reorder=reorder,
                            plot=plot(title=r"$\sigma_z$"),
                        ),
                    ],
                    "BSF-SPN": "BSF-SPOL",
                },
            ),
            NV(
                "I",
                "RAW_DATA",
                i(-1),
                reorder=reorder,
                plot=plot(negative=True, colormap=mymap, title="Total"),
            ),
        ],
        cls=BSFDefinition,
        name="BSF",
        info="BSF output file",
    )

    return definition


definition = create_definition()
