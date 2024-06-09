""" In this file the general plotting routines are present. """
import matplotlib
from matplotlib import rc_context
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, CenteredNorm, LogNorm, Normalize
from typing import Optional, Callable
import numpy as np
import functools
from matplotlib.colors import ListedColormap

from ..common.decorators import add_to_signature


def normalize_rc_params(params):
    out = {}
    for k,v in params.items():
        if isinstance(v, dict):
           v=normalize_rc_params(v)
           for kk,vv in v.items():
               out[k + '.' + kk] = vv
        else:
           out[k]=v
    return out


rc_params = normalize_rc_params({
    'font' : {'family' :'serif',
             'serif'  :['cms10'],
             'size'   : 16,
             'weight' :'normal'},
    'text' : { 'usetex' :True }
})


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


def create_rc_context(latex:bool=True):
    """
    Create the context that sets defaults for plotting
    """
    params = rc_params
    if not latex:
       params['text.usetex'] = False

    return rc_context(params)


def single_plot(fn:Callable, *args, filename:Optional[str]=None, show:Optional[bool]=None, dpi=600, latex=True, figsize=(6,4), callback=None, **kwargs):
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
      plt.subplots_adjust(left=0.15,right=0.95,bottom=0.17,top=0.93)
      fn(*args, **kwargs, axis=ax)
      if callback:
          callback(ax)
      finish_plot(filename, show, dpi)


def finish_plot(filename:Optional[str]=None, show:Optional[bool]=None, dpi=600):
     """
     Show the plot and/or save it to the given file
     """
     if show is None:
        show=filename is None
     if show:
        plt.show()
     if filename:
        plt.savefig(filename, dpi=dpi)


def auto_range(rng, data):
    """
    Fill the missing value in the given range by the data.

    >>> auto_range( (None, None), [2,5,-3,7] )
    (-3, 7)
    >>> auto_range( (None, 4), [2,5,-3,7] )
    (-3, 4)
    >>> auto_range( (2, 4), [2,5,-3,7] )
    (2, 4)
    """
    if rng is None:
       return ( np.min(data), np.max(data) )
    return (
           rng[0] if rng[0] is not None else np.min(data),
           rng[1] if rng[1] is not None else np.max(data),
    )


def plotting_function(func):
    """ Decorator, that 'completes' the given function that just draw into a
    matplolib axis.
    The completed function will have a few more arguments. One of them is
    ``axis``. If it is given, the plot is just drawn to the axis. If not,
    a plot is created, the function is called to draw into the plot, and
    then the plot is either showed or saved, according to the rest of the added
    arguments
    """

    @add_to_signature(func)
    @functools.wraps(func)
    def plot_function(*args, filename=None, show=None, dpi=600, latex=True, figsize=(6,4), callback=None, axis=None, **kwargs):
        if axis:
           func(*args, axis=axis, **kwargs)
           if callback:
             callback(axis)
        else:
           single_plot(func, filename=filename, show=show, dpi=dpi, latex=latex, callback=None, figsize=figsize, *args, **kwargs)
    return plot_function


def set_up_common_plot(axis, title=None, xlabel=None, ylabel=None, xticklabels=None, yticklabels=None, xticks=None, yticks=None):
   loc = locals()
   """
   This functions just set the properties of an matplotlib axis, that are common across various plots.
   """
   args = { n: loc[n]
            for n in ('xlabel', 'ylabel', 'xticks', 'yticks', 'xticklabels', 'yticklabels', 'title')
            if loc[n] is not None }
   for name in args:
       if args[name] is not None:
           getattr(axis, 'set_' + name)(args[name])


@plotting_function
@add_to_signature(set_up_common_plot, prepend=True, kwargs=True)
def colormesh(x,y,c, xrange=None, yrange=None, colormap=None, show_zero_line=False, axis=None, mode=False, norm=None, vmin=None, vmax=None, colorbar=False, **kwargs):
   """
   Plot 3D data by assigning colors to 2D grid. See matplotlib.pyplot.pcolormesh
   """
   set_up_common_plot(axis, **kwargs)
   if 'cmap' in kwargs:
      colormap = kwargs['cmap']

   if mode == 'centered':
       if norm == 'log':
           colormap = colormap or 'RdBu_r'
           norm = SymLogNorm(linthresh=1e-12,vmax=vmax)  # vmin=c.min(), vmax=c.max())
       elif norm =='lin':
           colormap = colormap or 'seismic'
           norm = CenteredNorm(vmax=vmax)
   else:
       colormap = colormap or 'BuPu'
       if norm == 'log':
           norm=LogNorm(vmin=1e-8, vmax=vmax)
       elif norm=='lin':
           norm=Normalize(vmin=0. if mode == 'from_zero' else None, vmax=vmax)

   axis.set_xlim(auto_range(xrange, x))
   axis.set_ylim(auto_range(yrange, y))
   axis.pcolormesh(x,y,c,cmap=colormap,shading='gouraud', norm=norm, vmin=None, vmax=None)
   if show_zero_line:
       opts = {
           'lw' : 1.,
           'color' : 'black'
       }
       if isinstance(show_zero_line, dict):
          opts.update(show_zero_line)
       elif show_zero_line is not True:
          opts['lw'] = show_zero_line
       axis.plot(axis.get_xlim(),[0,0],**opts)

   if colorbar:
     axis.figure.colorbar(plt.cm.ScalarMappable(matplotlib.colors.Normalize(vmin=vmin, vmax=vmax), cmap=colormap), ax=axis)


class Multiplot:
  """ This class can be used for plotting more plots into one resulting image/window. """

  def __init__(self, layout, figsize=(6,4), latex=True, updown_layout=False):
      self.fig, self.axes = plt.subplots(figsize=figsize, nrows=layout[0], ncols=layout[1])
      plt.subplots_adjust(left=0.12,right=0.95,bottom=0.17,top=0.90, hspace=0.75, wspace=0.5)
      self.free_axes = self.axes.ravel(order='F' if not updown_layout else 'C')
      self.free_axes = [ i for i in self.free_axes[::-1] ]

  def __iter__(self):
      axis = self.free_axes.pop()
      while axis:
          yield axis
          axis = self.free_axes.pop()

  def plot(self, option, plot_info=None, **kwargs):
      if not plot_info:
         plot_info = option._definition.plot
      plot_info(option, axis=self.free_axes.pop(), **kwargs)

  def finish(self, filename:Optional[str]=None, show:Optional[bool]=None, dpi=600):
      """
      Show the prepared plots or save them

      Parameters
      ----------
      filename
        The filename where the plot should be stored

      show
        If True, always show the plot.
        If False, never.
        If None, show it, if the filename is not set.

      dpi
        Dpi for generated plot
      """
      for i in self.free_axes:
          i.set_visible(False)
      finish_plot(filename, show, dpi)


def change_default_kwargs(f, **kwargs):
    """ Return the same function, with default kwargs changed """
    out = functools.partial(f, **kwargs)
    functools.update_wrapper(out, f)
    return out
