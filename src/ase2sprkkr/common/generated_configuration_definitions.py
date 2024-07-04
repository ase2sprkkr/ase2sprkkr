""" Definition of value in configuration/output files,
that is generated from other values """

from .configuration_definitions import RealItemDefinition
from .decorators import add_to_signature
from .options import Option
import copy
import numpy as np
from typing import Union


class BaseGeneratedValueDefinition(RealItemDefinition):
  """ Base class for all generated values. It just set
  that it is generated. """

  is_generated = True
  """ Generated value - the value is computed from other values """
  is_stored = False
  """ This property sets, that this Value/Option is generated. """

  result_class = Option
  """ The generated Values creates :class:`Option` """

  _grammar = None
  """ Do not generated grammar, since this item is not readed, but computed from other values. """

  def enrich(self, option):
      """ By default, generated values recieve no enhancement from
      its definition """
      pass

  def __repr__(self):
       return f"<{self.name} (generated)>"


class GeneratedValueDefinition(BaseGeneratedValueDefinition):

   @add_to_signature(RealItemDefinition.__init__, prepend=True)
   def __init__(self, name, getter, setter=None, **kwargs):
       super().__init__(name, **kwargs)
       self.getter = getter
       self._setter = setter

   @property
   def setter(self):
       if not self._setter:
           raise ValueError("Setting the value(s) of {self.name} is not allowed")
       return self._setter


class NumpyViewDefinition(BaseGeneratedValueDefinition):
   """
   Values described by this description are possibly reshaped views into a large
   "raw data array"

   Parameters
   ----------
   name
     Name of the resulting variable

   data
     The source variable, from which the data are taken

   selector
     The selector, which select the data to be viewed. Any slice
     or simple numpy index is allowed.

   shape
     The data will be reshaped to given shape. The dimension can be given
     either by number, or by names of other container variables. E.g.
     ``('NE', 5)``

   transpose
    Transpose the source data before returning or indexing.
    If reorder is given, this settings has no effect.

   reorder
    Reorder the axes after reshaping. Argument is the array
    of axes order, e.g. (2,0,1) shifts the last axis to be the first.

   transform_key
    Transform function for the keys idnexing the array
    (e.g. the string name can be transformed to a propper numerical index)

   plot
    PlotInfo object that defines how the results are plotted
   """

   def data_description(self, verbose:Union[bool,str]=False, show_hidden=False, prefix:str=''):
        if self.shape:
            shape=f"({self.shape})"
        else:
            shape=""

        out = f"{prefix}{self.name} : view of {self.data}{shape}"
        return out

   @add_to_signature(RealItemDefinition.__init__, prepend=True)
   def __init__(self, name, data, selector=slice(None),
                shape=None, transpose=False, reorder=None,
                transform_key=None,
                plot=None,
                *args, **kwargs):
       super().__init__(name, *args, **kwargs)
       self.selector = selector
       self.shape = shape
       self.data = data
       self.plot = plot
       self.transform_key=transform_key
       if reorder:
           self.reorder=reorder
       else:
           self.reorder=transpose

   def determine_shape(self, container):
       """ Return the shape of the resulting array, possibly computed using
       properties of the other values in the container"""
       def get(i):
           if isinstance(i, str):
              return container[i]()
           return i
       return tuple([get(i) for i in self.shape])

   def source(self, container):
       if callable(self.selector):
          out = self.selector(container[self.data](), container)
       else:
          out=container[self.data]()[self.selector]
       if self.shape:
          out.shape = self.determine_shape(container)
       if self.reorder:
          out = np.transpose(out, axes=self.reorder if self.reorder is not True else None)
       return out

   def getter(self, container, key=None):
       out=self.source(container)
       if key is not None:
          if self.transform_key:
             key=self.transform_key(key, container)
          out=out[key]
       return out

   def setter(self, container, value, key=slice(None)):
       if self.transform_key:
          key=self.transform_key(key, container)
       self.source(container)[key]=value

   def enrich(self, option):
       if self.plot:
         option.plot = lambda **kwargs: self.plot(option, **kwargs)
         option.plot.__doc__ = " Plot the data."

   def copy_value(self, value, all_values=False):
       return copy.copy(value)
