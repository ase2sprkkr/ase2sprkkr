""" In this file the general plotting routines are present. """

from matplotlib import rc_context
import matplotlib.pyplot as plt
from typing import Optional, Callable
import numpy as np
import functools
from ..common.decorators import add_to_signature

def normalize_rc_params(params):
    out = {}
    for k,v in params.items():
        if isinstance(v, dict):
           v=normalize_rc_params(v)
           for kk,vv in v.items():
               out[k+'.'+kk] = vv
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


from matplotlib.colors import ListedColormap
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

def single_plot(fn:Callable, *args, filename:Optional[str]=None, show:Optional[bool]=None, dpi=600, latex=True, figsize=(6,4), **kwargs):
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
      finish_plot(filename, show, dpi)

def finish_plot(filename:Optional[str]=None, show:Optional[bool]=None, dpi=600):
     if show is None:
        show=filename is None
     if show:
        plt.show()
     if filename:
        plt.savefig(filename, dpi=dpi)


def auto_range(rng, data):
    if rng is None:
       return ( np.min(data), np.max(data) )
    return (
        data[0] if data[0] is not None else np.min(data),
        data[1] if data[1] is not None else np.max(data),
        )
def plotting_function(func):
    @add_to_signature(func)
    @functools.wraps(func)
    def plot_function(*args, filename=None, show=None, dpi=600, latex=True, axis=None, **kwargs):
        if axis:
           func(*args, axis=axis, **kwargs)
        else:
           single_plot(func, filename=filename, show=show, dpi=dpi, latex=latex, *args, **kwargs)
    return plot_function

def common_plot(title=None, xlabel=None, ylabel=None, xticklabels=None, yticklabels=None, axis=None):
   l = locals()
   args = { n: l[n]
            for n in ('xlabel', 'ylabel', 'xticklabels', 'yticklabels', 'title')
            if l[n] is not None }
   for name in args:
       if args[name] is not None:
           getattr(axis, 'set_'+name)(args[name])


@plotting_function
@add_to_signature(common_plot, prepend=True)
def colormesh(x,y,c, xrange=None, yrange=None, colormap=None, show_zero_line=False, axis=None, **kwargs):
   """
   Plot 3D data by assigning colors to 2D grid. See matplotlib.pyplot.pcolormesh
   """
   common_plot(**kwargs, axis=axis)
   colormap = colormap or 'Blues'
   axis.set_xlim(auto_range(xrange, x))
   axis.set_ylim(auto_range(yrange, y))
   axis.pcolormesh(x,y,c,cmap=colormap,shading='gouraud',
                 vmin=c.min(), vmax=c.max())
   if show_zero_line:
       axis.plot(axis.get_xlim(),[0,0],color='black',lw=1)


class Multiplot:

  def __init__(self, layout, figsize=(6,4), latex=True, updown_layout=False):
      self.fig, self.axes = plt.subplots(figsize=(6,4), nrows=layout[0], ncols=layout[1])
      plt.subplots_adjust(left=0.12,right=0.95,bottom=0.17,top=0.90, hspace=0.75, wspace=0.5)
      self.free_axes = self.axes.ravel(order='F' if not updown_layout else 'C')
      self.free_axes = [ i for i in self.free_axes[::-1] ]


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
      finish_plot(filename, show, dpi)

class PlotInfo:
  """
  Object that describe, how the values are plotted (the method and default values
  for plotting)
  """

  def __init__(self, axes, method=None, **kwargs):
      self.kwargs = kwargs
      self.axes = axes
      self.method = method

  def __call__(self, option, **kwargs):

      values = option()
      tmp = self.kwargs.copy()
      tmp.update(kwargs)
      kwargs = tmp
      if 'args' in kwargs and option.name in kwargs['args']:
         kwargs.update(kwargs['args'][option.name])

      axes = kwargs.get('axes', self.axes)
      method = kwargs.get('method', self.method)
      if method is None:
         if values.ndim == 2: method = colormesh

      kwargs = { i:j for i,j in kwargs.items() if i not in ('axes', 'method', 'args') }

      xdata = []
      for i,axe in zip(axes, ['x', 'y', 'z']):
          pname = axe + 'label'
          if isinstance(i, tuple):
              i, label = i
          else:
              label = i
          if label and not pname in kwargs:
              kwargs[pname] = label
          if isinstance(i, str):
             i = option._container[i]()
          xdata.append(i)

      if 'title' in kwargs and kwargs['title'] == False:
         del kwargs['title']
      elif 'title' not in kwargs:
         kwargs['title'] = option.info if option.info else option.name

      method(*xdata, values, **kwargs)
