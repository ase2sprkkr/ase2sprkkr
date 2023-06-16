""" In this file the general plotting routines are present. """

from matplotlib import rc_context
import matplotlib.pyplot as plt
from typing import Optional, Callable
import numpy as np

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


def do_plotting(fn:Callable, filename:Optional[str]=None, show:Optional[bool]=None, dpi=600, latex=True):
    """
    Do the actual plotting (sets the environment, save the plot etc.)

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
    params = rc_params
    if not latex:
       params['text.usetex'] = False

    with rc_context(params):
       fn()

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

def colormesh(x,y,c, filename=None, show=None, xrange=None, yrange=None, colormap=None, dpi=600,
              xlabel=None, ylabel=None, xticklabels=None, yticklabels=None,
              show_zero_line=False, latex=True
              ):
   """
   Plot 3D data by assigning colors to 2D grid. See matplotlib.pyplot.pcolormesh
   """
   l = locals()
   args = { n: l[n]
            for n in ('xlabel', 'ylabel', 'xticklabels', 'yticklabels')
            if l[n] is not None }
   colormap = colormap or 'Blues'

   def plot():
       fig, ax = plt.subplots(figsize=(6,4))
       plt.subplots_adjust(left=0.15,right=0.95,bottom=0.17,top=0.93)
       ax.set_xlim(auto_range(xrange, x))
       ax.set_ylim(auto_range(yrange, y))
       ax.pcolormesh(x,y,c,cmap=colormap,shading='gouraud',
                     vmin=c.min(), vmax=c.max())
       for name in args:
           getattr(ax, 'set_'+name)(args[name])
       if show_zero_line:
           ax.plot(xrange,[0,0],color='black',lw=1)


   do_plotting(plot, filename, show, dpi, latex)


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
      axes = kwargs.get('axes', self.axes)
      method = kwargs.get('method', self.method)
      if method is None:
         if values.ndim == 2: method = colormesh

      kwargs = { i:j for i,j in kwargs.items() if i not in set('axes', 'method') }

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

      method(*xdata, values, **kwargs)
