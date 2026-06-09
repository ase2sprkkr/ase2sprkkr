"""In this file the general plotting routines are present."""

import matplotlib
from matplotlib import rc_context
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, CenteredNorm, LogNorm, Normalize
from matplotlib.colors import ListedColormap

from typing import Optional
import numpy as np
import os
import numbers
import re
import functools
from contextlib import contextmanager

from ..common.decorators import add_to_signature


def normalize_rc_params(params):
    out = {}
    for k, v in params.items():
        if isinstance(v, dict):
            v = normalize_rc_params(v)
            for kk, vv in v.items():
                out[k + "." + kk] = vv
        else:
            out[k] = v
    return out


rc_params = normalize_rc_params(
    {"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10, "font.weight": "normal"}
)


def combine_colormaps(cmap1, cmap2, n1, n2, index1, index2):

    # Get the selected colors from each colormap
    colors1 = cmap1(np.linspace(index1[0], index1[1], n1))
    colors2 = cmap2(np.linspace(index2[0], index2[1], n2))

    # Combine the selected colors from the two colormaps
    combined_colors = np.vstack((colors1, colors2))

    # Create the new colormap
    combined_cmap = ListedColormap(combined_colors)

    return combined_cmap


def combined_colormap(range1=(0.5, 0), range2=(0.15, 1), n1=8000, n2=15000, cmap1=plt.cm.bwr, cmap2=plt.cm.jet):
    # Create the new colormap
    return combine_colormaps(cmap1, cmap2, n1, n2, range1, range2)


def create_rc_context(latex: Optional[bool] = None):
    """
    Create the context that sets defaults for plotting
    """
    params = rc_params
    if latex is not None:
        params["text.usetex"] = latex
    return rc_context(params)


@contextmanager
def single_plot(
    filename: Optional[str] = None, show: Optional[bool] = None, window_title=None, dpi=600, latex=None, figsize=(6, 4), **kwargs
):
    """
    Creates single plot according to the given function a either show it or save it.

    Parameters
    ----------
    fn
      Function that sets the plot

    filename
      The filename where the plot should be stored

    show
      If True, always show the plot.
      If False, never.
      If None, show it, if the filename is not set.
    """
    with create_rc_context(latex):
        fig, ax = plt.subplots(figsize=figsize)
        if window_title:
            fig.canvas.manager.set_window_title(window_title)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.17, top=0.93)
        yield ax, kwargs
        finish_plot(fig, filename, show, dpi)


def finish_plot(fig, filename: Optional[str] = None, show: Optional[bool] = None, dpi=600, layout=False):
    """
    Show the plot and/or save it to the given file
    """
    if layout == "tight":
        plt.tight_layout()
    if show is None:
        show = filename is None
    if filename:
        plt.savefig(filename, dpi=dpi)

    if layout == "constrained":
        try:
            fig.get_layout_engine().set(
                w_pad=0.1,  # horizontal outer padding
                h_pad=0.1,  # vertical outer padding
                wspace=0.05,
                hspace=0.05,
            )
        except AttributeError:
            fig.set_constrained_layout_pads(
                w_pad=0.1,  # horizontal outer padding
                h_pad=0.1,  # vertical outer padding
                wspace=0.05,
                hspace=0.05,
            )
    if show:
        plt.show()


def auto_range(rng, data, eps=1e-4):
    """
    Fill the missing value in the given range by the data.

    >>> if np.__version__ > "2.0":
    ...     np.set_printoptions(legacy="1.25")
    >>> auto_range((None, None), [2, 5, -3, 7])
    (-3, 7)
    >>> auto_range((None, 4), [2, 5, -3, 7])
    (-3, 4)
    >>> auto_range((2, 4), [2, 5, -3, 7])
    (2, 4)
    """
    if rng is None:
        out = (np.min(data), np.max(data))
    else:
        out = (rng[0] if rng[0] is not None else np.min(data), rng[1] if rng[1] is not None else np.max(data))
    if out[0] == out[1]:
        val = out[0]
        if val == 0:
            out = (val - eps, val + eps)
        elif out[0] > 0:
            out = (max(val - eps, 0), val + eps)
        else:
            out = (val - eps, min(val + eps, 0))
    return out


def plotting_function(func):
    """Decorator, that 'completes' the given function that just draw into a
    matplolib axis.
    The completed function will have a few more arguments. One of them is
    ``axis``. If it is given, the plot is just drawn to the axis. If not,
    a plot is created, the function is called to draw into the plot, and
    then the plot is either showed or saved, according to the rest of the added
    arguments
    """

    @add_to_signature(func)
    @functools.wraps(func)
    def plot_function(
        *args, filename=None, show=None, dpi=600, latex=None, figsize=(6, 4), callback=None, axis=None, **kwargs
    ):
        if axis:
            func(*args, axis=axis, **kwargs)
            if callback:
                callback(axis)
        else:
            with single_plot(filename=filename, show=show, dpi=dpi, latex=latex, figsize=figsize) as (axis, _):
                func(*args, axis=axis, **kwargs)
                if callback:
                    callback(axis)

    return plot_function


def set_up_common_plot(
    axis, title=None, xlabel=None, ylabel=None, xticklabels=None, yticklabels=None, xticks=None, yticks=None, **kwargs
):
    loc = locals()
    """
   This functions just set the properties of an matplotlib axis, that are common across various plots.
   """
    args = {
        n: loc[n]
        for n in ("xlabel", "ylabel", "xticks", "yticks", "xticklabels", "yticklabels", "title")
        if n != "kwargs" and loc[n] is not None
    }
    kwargs.update(args)
    for name in kwargs:
        if not hasattr(axis, "set_" + name):
            raise ValueError(f"Axis has not set_{name} method, thus I don't know what to do with {name} argument")
        getattr(axis, "set_" + name)(kwargs[name])


@plotting_function
@add_to_signature(set_up_common_plot, prepend=True, kwargs=True)
def colormesh(
    x,
    y,
    c,
    xrange=None,
    yrange=None,
    colormap=None,
    show_zero_line=False,
    axis=None,
    mode=False,
    norm=None,
    vmin=None,
    vmax=None,
    colorbar=False,
    **kwargs,
):
    """
    Plot 3D data by assigning colors to 2D grid. See matplotlib.pyplot.pcolormesh
    """
    set_up_common_plot(axis, **kwargs)
    if "cmap" in kwargs:
        colormap = kwargs["cmap"]

    if mode == "centered":
        if norm == "log":
            colormap = colormap or "RdBu_r"
            norm = SymLogNorm(linthresh=1e-12, vmax=vmax)  # vmin=c.min(), vmax=c.max())
            vmax = None
        elif norm == "lin":
            colormap = colormap or "seismic"
            norm = CenteredNorm(vmax=vmax)
            vmax = None
    else:
        colormap = colormap or "BuPu"
        if norm == "log":
            if vmin is None:
                if vmax is not None and vmax > 0.0:
                    vmin = min(1e-8, vmax)
                else:
                    vmin = 1e-8
            if vmax is not None and vmax < vmin:
                vmax = vmin
            norm = LogNorm(vmin=vmin, vmax=vmax)
            vmax = None
            vmin = None
        elif norm == "lin":
            norm = Normalize(vmin=0.0 if mode == "from_zero" else None, vmax=vmax)
            vmax = None
            vmin = None

    axis.set_xlim(auto_range(xrange, x))
    axis.set_ylim(auto_range(yrange, y))
    axis.pcolormesh(x, y, c, cmap=colormap, shading="gouraud", norm=norm, vmin=vmin, vmax=vmax)
    if show_zero_line:
        opts = {"lw": 1.0, "color": "black"}
        if isinstance(show_zero_line, dict):
            opts.update(show_zero_line)
        elif show_zero_line is not True:
            opts["lw"] = show_zero_line
        axis.plot(axis.get_xlim(), [0, 0], **opts)

    if colorbar:
        axis.figure.colorbar(
            plt.cm.ScalarMappable(matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap=colormap), ax=axis
        )


class Multiplot:
    """This class can be used for plotting more plots into one resulting image/window."""

    def __init__(
        self,
        layout,
        figsize=(6, 4),
        latex=None,
        updown_layout=False,
        filename: Optional[str] = None,
        show: Optional[bool] = None,
        dpi=600,
        separate_plots=False,
        layout_kind="constrained",
        number_of_plots=None,
        **kwargs,
    ):
        self.separate_plots = separate_plots
        self.filename = filename
        self.show = show
        self.dpi = dpi
        self.latex = latex
        self.layout_kind = layout_kind  #'constrained', 'tight', 'adjust', {dict for adjust}

        if isinstance(layout, numbers.Integral):
            layout = min(layout, number_of_plots)
            layout = ((number_of_plots - 1) // layout + 1, layout)
        elif layout[0] is None:
            layout = ((number_of_plots - 1)  // layout[1] + 1, layout[1])
        elif layout[1] is None:
            layout = (layout[0], (number_of_plots -1) // layout[0] + 1)

        if callable(figsize):
            figsize = figsize(layout)

        if separate_plots:
            self.figsize = figsize
            self.number = layout[0] * layout[1]
            self.figure = None
        else:
            self.figure, self.axes = plt.subplots(
                figsize=figsize, nrows=layout[0], ncols=layout[1], constrained_layout=layout_kind == "constrained"
            )
            if layout_kind == "adjust" or isinstance(layout_kind, dict):
                adj = {"left": 0.12, "right": 0.95, "bottom": 0.17, "top": 0.90, "hspace": 0.75, "wspace": 0.4}
                if isinstance(layout_kind, dict):
                    adj.update(self.layout_kind)
                plt.subplots_adjust(**adj)

            self.free_axes = np.atleast_1d(self.axes).ravel(order="C" if not updown_layout else "F")
            self.free_axes = [i for i in self.free_axes[::-1]]

        self.specific_kwargs = {k: v for k, v in kwargs.items() if str(k).isnumeric()}
        for i in self.specific_kwargs:
            del kwargs[i]
        self.kwargs = kwargs
        self.specific_kwargs = {int(k): v for k, v in self.specific_kwargs.items()}
        self.index = 0

    def __enter__(self):
        if not self.separate_plots:
            self.context = create_rc_context(latex=self.latex)
            self.context.__enter__()
        return self

    def __exit__(self, type, value, traceback):
        if not self.separate_plots:
            for i in self.free_axes:
                i.set_visible(False)
            finish_plot(self.figure, self.filename, self.show, self.dpi, layout=self.layout_kind)
            return self.context.__exit__(type, value, traceback)
        else:
            if self.separate_plots != "show":
                show = self.show
                if show is None:
                    show = not self.filename
                if show:
                    plt.show()

    def plot(self, option, name=None, filename=None, plot_function=None, **kwargs):
        if not plot_function:
            if option.number_of_plots() > 1:
                 option.plot(filename=filename, multiplot=self, **kwargs)
                 return

            def plot_function(**kwargs):
                return option.plot(**kwargs)

        name = name or getattr(option, "name", None) or str(option)
        with self.new_axis(name=name, filename=filename, **kwargs) as (axis, kw):
            plot_function(axis=axis, **kw)

    @contextmanager
    def new_axis(self, name=None, filename=None, **kwargs):
        for k, v in self.kwargs.items():
             kwargs.setdefault(k, v)
        for k, v in self.specific_kwargs.get(self.index, {}).items():
             kwargs.setdefault(k, v)

        if self.separate_plots:
            if self.index >= self.number:
                raise StopIteration()

            def append_before_ext(filename: str, suffix: str) -> str:
                root, ext = os.path.splitext(filename)
                if ext:
                    return f"{root}{suffix}{ext}"
                else:
                    return f"{root}{suffix}"

            if filename is None:
                filename = self.filename
                if filename:
                    if name is None:
                        fname = str(self.index + 1)
                    else:
                        fname = re.sub(r'[<>:"/\\|?*\x00-\x1f\x7f ]', "_", name)
                        fname = re.sub(r"_+", "_", fname)

                    if "{name}" in filename:
                        filename = filename.replace("{name}", fname)
                    else:
                        filename = append_before_ext(filename, "_" + fname)

            with single_plot(
                filename,
                self.show if self.separate_plots == "each" else False,
                name,
                self.dpi,
                self.latex,
                self.figsize,
            ) as (axis, _):
                yield axis, kwargs
        else:
            try:
                axis = self.free_axes.pop()
            except IndexError:
                raise StopIteration()
            yield axis, kwargs
        self.index += 1

    def __iter__(self):
        while True:
            with self.new_axis() as out:
                yield out


def change_default_kwargs(f, **kwargs):
    """Return the same function, with default kwargs changed"""
    out = functools.partial(f, **kwargs)
    functools.update_wrapper(out, f)
    return out
