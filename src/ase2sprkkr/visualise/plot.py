""" In this file the general plotting routines are present. """
import matplotlib
from matplotlib import rc_context
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, CenteredNorm, LogNorm, Normalize
from matplotlib.colors import ListedColormap

from typing import Optional, Callable
import numpy as np
import os
import re
import functools
from contextlib import contextmanager

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
    "font.family" : "serif",
    "mathtext.fontset": "cm",
    "font.size": 10,
    "font.weight": "normal"
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


def create_rc_context(latex:Optional[bool]=None):
    """
    Create the context that sets defaults for plotting
    """
    params = rc_params
    if latex is not None:
        params['text.usetex'] = latex
    return rc_context(params)

@contextmanager
def single_plot(filename:Optional[str]=None, show:Optional[bool]=None, dpi=600, latex=None, figsize=(6,4)):
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
      yield ax
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

    >>> if np.__version__ > '2.0': np.set_printoptions(legacy='1.25')
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
    def plot_function(*args, filename=None, show=None, dpi=600, latex=None, figsize=(6,4), callback=None, axis=None, **kwargs):
        if axis:
           func(*args, axis=axis, **kwargs)
           if callback:
             callback(axis)
        else:
           with single_plot(func, filename=filename, show=show, dpi=dpi, latex=latex, figsize=figsize) as axis:
                func(*args, axis=axis, **kwargs)
                if callback:
                    callback(axis)
    return plot_function


def set_up_common_plot(axis, title=None, xlabel=None, ylabel=None, xticklabels=None, yticklabels=None, xticks=None, yticks=None, **kwargs):
   loc = locals()
   """
   This functions just set the properties of an matplotlib axis, that are common across various plots.
   """
   args = { n: loc[n]
            for n in ('xlabel', 'ylabel', 'xticks', 'yticks', 'xticklabels', 'yticklabels', 'title')
            if n != 'kwargs' and loc[n] is not None }
   kwargs.update(args)
   for name in kwargs:
       if not hasattr(axis, 'set_' + name):
           raise ValueError(f"Axis has not set_{name} method, thus I don't know what to do with {name} argument")
       getattr(axis, 'set_' + name)(kwargs[name])


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
           vmax = None
       elif norm =='lin':
           colormap = colormap or 'seismic'
           norm = CenteredNorm(vmax=vmax)
           vmax = None
   else:
       colormap = colormap or 'BuPu'
       if norm == 'log':
           norm=LogNorm(vmin=1e-8, vmax=vmax)
           vmax = None
       elif norm=='lin':
           norm=Normalize(vmin=0. if mode == 'from_zero' else None, vmax=vmax)
           vmax = None

   axis.set_xlim(auto_range(xrange, x))
   axis.set_ylim(auto_range(yrange, y))
   axis.pcolormesh(x,y,c,cmap=colormap,shading='gouraud', norm=norm, vmin=vmin, vmax=vmax)
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

  def __init__(self, layout, figsize=(6,4), latex=None, updown_layout=False,
               filename:Optional[str]=None, show:Optional[bool]=None, dpi=600,
               separate_plots=False, adjust={},
               **kwargs):
      self.separate_plots = separate_plots
      self.filename = filename
      self.show = show
      self.dpi = dpi
      self.latex = latex

      if separate_plots:
          self.figsize = figsize
          self.number = layout[0] * layout[1]
          self.figure = None
      else:
          self.figure, self.axes = plt.subplots(figsize=figsize, nrows=layout[0], ncols=layout[1])
          adj = {'left': 0.12, 'right': 0.95, 'bottom': 0.17, 'top': 0.90, 'hspace': 0.75, 'wspace': 0.5}
          adj.update(adjust)
          plt.subplots_adjust(**adjust)
          self.free_axes = self.axes.ravel(order='F' if not updown_layout else 'C')
          self.free_axes = [ i for i in self.free_axes[::-1] ]

      self.specific_kwargs = { k:v for k,v in kwargs.items() if str(k).isnumeric() }
      for i in self.specific_kwargs:
          del kwargs[i]
      self.kwargs = kwargs
      self.specific_kwargs = { int(k):v for k,v in self.specific_kwargs.items() }
      self.index = 0

  def __enter__(self):
      if not self.separate_plots:
          self.context = create_rc_context(latex=self.latex)
          self.context.__enter__()
      return self

  def __exit__(self,type, value, traceback):
      if not self.separate_plots:
          for i in self.free_axes:
              i.set_visible(False)
          finish_plot(self.filename, self.show, self.dpi)
          return self.context.__exit__(type, value, traceback)
      else:
          if self.separate_plots!='show':
              show = self.show
              if show is None:
                  show = not self.filename
              if show:
                  plt.show()

  def plot(self, option, name=None, plot_function=None, **kwargs):
      name = name or getattr(option, 'name', None) or str(option)
      with self.new_axis(name=name) as axis:
          kw = self.kwargs.copy()
          kw.update(self.specific_kwargs.get(self.index, {}))
          kw.update(kwargs)
          if not plot_function:
              plot_function = lambda **kwargs: option.plot(**kwargs)
          plot_function(axis=axis, **kw)

  @contextmanager
  def new_axis(self, name=None):
      if self.separate_plots:
          if self.index >= self.number:
             raise StopIteration()

          def append_before_ext(filename: str, suffix: str) -> str:
              root, ext = os.path.splitext(filename)
              if ext:
                  return f"{root}{suffix}{ext}"
              else:
                  return f"{root}{suffix}"

          filename = self.filename
          if filename:
              if name is None:
                  fname = str(self.index+1)
              else:
                  fname = re.sub(r'[<>:"/\\|?*\x00-\x1f\x7f ]', '_', name)
                  fname = re.sub(r"_+", '_', fname)

              filename = append_before_ext(filename, '_' + fname)

          with single_plot(filename, self.show if self.separate_plots=='each' else False,
                           self.dpi, self.latex, self.figsize) as axis:
               if name:
                   axis.figure.canvas.manager.set_window_title(name)
               yield axis
      else:
          try:
              yield self.free_axes.pop()
          except IndexError:
              raise StopIteration()
      self.index += 1

  def __iter__(self):
      while True:
          with self.new_axis() as ax:
              yield ax



def change_default_kwargs(f, **kwargs):
    """ Return the same function, with default kwargs changed """
    out = functools.partial(f, **kwargs)
    functools.update_wrapper(out, f)
    return out
