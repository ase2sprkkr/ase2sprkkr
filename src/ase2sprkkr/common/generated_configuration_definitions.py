""" Definition of value in configuration/output files,
that is generated from other values """

from .configuration_definitions import BaseDefinition
from .decorators import add_to_signature
from functools import partial
from .options import Option
import copy

class BaseGeneratedValueDefinition(BaseDefinition):
  """ Base class for all generated values. It just set
  that it is generated. """

  is_generated = True
  """ This property sets, that this Value/Option is generated. """

  result_class = Option

  def enrich(self, option):
      """ By default, generated values recieve no enhancement from
      its definition """
      pass

class GeneratedValueDefinition(BaseGeneratedValueDefinition):

   @add_to_signature(BaseDefinition.__init__, prepend=True)
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
   """

   @add_to_signature(BaseDefinition.__init__, prepend=True)
   def __init__(self, name, data, selector=slice(None),
                shape=None, transpose=False, transform_key=None,
                plot=None,
                *args, **kwargs):
       super().__init__(name, *args, **kwargs)
       self.selector = selector
       self.shape = shape
       self.data = data
       self.plot = plot
       self.transform_key=transform_key
       self.transpose = transpose

   def determine_shape(self, container):
       """ Return the shape of the resulting array, possibly computed using
       properties of the other values in the container"""
       def get(i):
           if isinstance(i, str):
              return container[i]()
           return i
       return tuple([get(i) for i in self.shape])

   def source(self, container):
       out=container[self.data]()[self.selector]
       if self.shape:
          out.shape = self.determine_shape(container)
       if self.transpose:
          out = out.T
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
