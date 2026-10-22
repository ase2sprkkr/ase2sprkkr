"""Structured four-column KKRSPEC BAND results and band-structure plots.

BAND rows are discrete accepted Bloch roots rather than regularly gridded
intensities. The native first coordinate is kz for NT=1 and k_parallel for
NT>1; the second coordinate is energy relative to the Fermi level; the final two
columns contain the real and imaginary parts of complex T_EIG.
"""

from typing import Optional

from ..result_files_definitions import (
    ResultFileDefinition,
    ResultFileValueDefinition as V,
    create_result_file_definition,
)
from ..result_files import CommonResultFile
from ...common.grammar_types import NumpyArray, Prefixed
from ...common.generated_configuration_definitions import NumpyViewDefinition as NV
from ...gui.plot import single_plot


class BandSPCResultFile(CommonResultFile):
    """Native BAND roots, exposed as a single energy-versus-wavevector plot."""

    def T_EIG(self):
        """Complex transmission eigenvalues written as two real columns."""
        return self.T_EIG_REAL() + 1j * self.T_EIG_IMAG()

    def plot(
        self,
        *,
        filename: Optional[str] = None,
        show: Optional[bool] = None,
        axis=None,
        figsize=(7, 5),
        dpi=600,
        latex=None,
        show_fermi=True,
        xlabel=None,
        ylabel=r"$E-E_{\mathrm{F}}$ (eV)",
        **kwargs,
    ):
        """Plot accepted band roots as points, never join unrelated eigenstates.

        For NT=1 the native K_EIG coordinate represents real kz; otherwise
        it represents k_parallel. The complex transmission eigenvalue is not used for plotting.
        """
        data = self.RAW_DATA()
        if xlabel is None:
            xlabel = (
                r"$k_z$ (native units)" if self.NT() == 1
                else r"$k_{\parallel}$ (native units)"
            )

        style = {"s": 6, "marker": ".", "linewidths": 0}
        style.update(kwargs)
        with single_plot(
            filename=filename,
            show=show,
            axis=axis,
            figsize=figsize,
            dpi=dpi,
            latex=latex,
        ) as (ax, _):
            ax.scatter(data[:, 0], data[:, 1], **style)
            if show_fermi:
                ax.axhline(0, color="0.5", linestyle="--", linewidth=0.8)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        return ax


class BandSPCDefinition(ResultFileDefinition):
    result_class = BandSPCResultFile

    @property
    def extension(self):
        return "spc"


definition = create_result_file_definition(
    "BAND",
    [
        V("NT", int),
        V("NP", int),
        V("COMMENT", Prefixed("#"), name_in_grammar=False),
        V(
            "RAW_DATA",
            NumpyArray(
                shape=(-1, 4),
                written_shape=(-1, 4),
                trivial_parsing=float
            ),
            name_in_grammar=False,
        ),
        NV("K_EIG", "RAW_DATA", (slice(None), 0)),
        NV("E_EIG", "RAW_DATA", (slice(None), 1)),
        NV("T_EIG_REAL", "RAW_DATA", (slice(None), 2)),
        NV("T_EIG_IMAG", "RAW_DATA", (slice(None), 3)),
    ],
    name="BAND_SPC",
    cls=BandSPCDefinition,
    info="KKRSPEC BAND output: native wavevector coordinate, energy relative to EF, and the real and imaginary components of T_EIG.",
)
