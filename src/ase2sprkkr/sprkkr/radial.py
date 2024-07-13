import numpy as np
import copy

from ..common.decorators import cached_property


class RadialValue:

  def __init__(self, value, mesh):
      self._value = np.asarray(value)
      self._mesh = mesh

  @cached_property
  def interpolator(self):

      def inter(x):
          return self._mesh.interpolator(x)
      ndim = self._value.ndim
      if ndim == 1:
          return inter(self._value)

      shape = self._value.shape[0]

      inters = [ inter(self._value[i]) for i in range(shape) ]

      def multiinter(x):
          if isinstance(x, np.ndarray):
              sx = (shape, ) + x.shape
          else:
              sx = (shape, )
          out = np.empty(sx)
          for i in range(shape):
              out[i] = inters[i](x)
          return out

      return multiinter

  @property
  def mesh(self):
      return self._mesh

  @mesh.setter
  def mesh(self, mesh):
      if mesh is self._mesh:
          return
      value = self.interpolate(mesh.coors)
      self._mesh=mesh
      self._value = value

  def copy(self, mesh=None):
      out = copy.copy(self)
      if mesh is None or self._mesh is mesh:
          out._value = np.copy(self._value)
      else:
          out.mesh = mesh
      return out

  def for_mesh(self, mesh):
      if mesh is self.mesh:
          return self
      return self.copy(mesh)

  def interpolate(self, coors):
      return self.interpolator(coors)

  def __call__(self, at):
      return self.interpolate(at)

  @property
  def raw_value(self):
      return self._value


class RadialCharge(RadialValue):
      pass


class RadialPotential(RadialValue):

  @staticmethod
  def from_munchen(value, mesh, z):
      return RadialPotential(value, mesh, z)

  @staticmethod
  def from_julich(value, mesh, z):
      v = np.empty_like(value)
      out = RadialPotential(v, mesh, z)
      out.julich = value
      return out

  def __init__(self, value, mesh, z):
      super().__init__(value, mesh)
      self.z = z

  @property
  def vt(self):
      return self._value[0]

  @vt.setter
  def vt(self, val):
      self._value[0,:]=val
      self._clear()

  @property
  def bt(self):
      return self._value[1]

  @bt.setter
  def bt(self, val):
      self._value[1,:]=val
      self._clear()

  @property
  def munchen(self):
      """ Potential in Munchen format """
      return self._value

  @munchen.setter
  def munchen(self, value):
      self._value = value
      self._clear()

  @property
  def julich(self):
      """ Potential in Julich format """
      out=np.empy_like(self._value)
      out[0]=self.up
      out[1]=self.down

  @julich.setter
  def julich(self, value):
      self._value[0] = 0.5 * (value[0] + value[1]) - 2 * self.zr
      self._value[1] = 0.5 * (value[0] - value[1])
      self._clear()

  @property
  def up(self):
      return self._value[0] + self._value[1] + 2 * self.zr

  @property
  def down(self):
      return self._value[0] - self._value[1] + 2 * self.zr

  def _clear(self):
      """ Clear cache on a value change """
      if 'interpolator' in self.__dict__:
          del self.__dict__['interpolator']

  @cached_property
  def zr(self):
      return self._z / self._mesh.coors
