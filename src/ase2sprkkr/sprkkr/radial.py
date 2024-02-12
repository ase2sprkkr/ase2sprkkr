import numpy as np
import copy

from ..common.decorators import cached_property


class RadialValue:

  def __init__(self, value, mesh):
      self._value = value
      self._mesh = mesh

  def recompute_for_mesh(self, mesh):
      out = copy.copy(self)
      out._mesh = mesh
      out._value = np.empty_like(out._value)
      out._value[0] = np.interp(mesh.coors, self._mesh.coors, self._value[0])
      out._value[1] = np.interp(mesh.coors, self._mesh.coors, self._value[1])


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

  @property
  def bt(self):
      return self._value[1]

  @bt.setter
  def bt(self, val):
      self._value[1,:]=val

  @property
  def munchen(self):
      """ Potential in Munchen format """
      return self._value

  @munchen.setter
  def munchen(self, value):
      self._value = value

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

  @property
  def up(self):
      return self._value[0] + self._value[1] + 2 * self.zr

  @property
  def down(self):
      return self._value[0] - self._value[1] + 2 * self.zr

  @cached_property
  def zr(self):
      return self._z / self._mesh.coors
