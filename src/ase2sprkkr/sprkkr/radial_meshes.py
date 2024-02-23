""" Definitions of radial meshes used by SPR-KKR. """
import copy
import numpy as np
import math
from scipy.interpolate import CubicSpline


class Mesh:
  """
  A base class for SPR-KKR radial meshes.
  """
  @staticmethod
  def default():
      """ SPR-KKR computes the mesh itself, if zeros are given """
      return ExponentialMesh(1e-6,2e-2,0,0.,721, 0.)

  def interpolate(self, coors, value):
      """ Interpolate values in the given coors to the mesh coors.
      """
      return self.interpolator(self.coors, value)(coors)

  def interpolator(self, values):
      return self.interpolator_for_coors(self.coors, values)

  @staticmethod
  def interpolator_for_coors(coors, values):
      v1 = abs(values[0])
      v2 = abs(values[-1])

      if v1 == 0.:
          r_method = False
      elif v2 == 0.:
          r_method = True
      else:
          c1 = coors[0]
          c2 = coors[-1]
          d = (v1 - v2)
          dr = (v1 * c1 - v2 * c2)
          r_method = abs(d) > abs(dr)

      if r_method:
          spl = CubicSpline(coors, values * coors, extrapolate=True)
          return lambda x: spl(x) / x

      return CubicSpline(coors, values)


def _clearing_property(name):

    attr = '_' + name

    def reader(self):
        return getattr(self, attr)

    def writer(self, value):
        setattr(self, attr, value)
        self._clear()

    out = property(reader)
    out = out.setter(writer)
    return out


class ExponentialMesh(Mesh):
  """
  Radial mesh definition for an atomic site.
  """

  def __init__(self, r1, dx, jrmt, rmt, jrws, rws):
      """
      Parameters
      ----------
      r1 : float
        First mesh item

      dx : float
        Multiplier: r_n = r_{n-1} * dx

      jrmt: int

      rmt: float

      jrws: int

      rws: float
      """
      self._r1 = r1
      self._dx = dx
      self._jrmt = jrmt
      self._rmt = rmt
      self._jrws = jrws
      self._rws = rws
      self._coors = None

  r1 = _clearing_property('r1')
  dx = _clearing_property('dx')
  jrmt = _clearing_property('jrmt')
  rmt = _clearing_property('rmt')
  jrws = _clearing_property('jrws')
  rws = _clearing_property('rws')

  def __len__(self):
      return self.jrws

  def __bool__(self):
      return True

  def __getitem__(self, i):
      return self.coors[i]

  def _clear(self):
      self._coors = None

  @property
  def coors(self):
      if self._coors is None:
          self._coors = np.empty(self.jrws)
          self._coors[0] = r = self.r1
          dx = math.exp(self.dx)
          for i in range(1, self.jrws):
              self._coors[i] = r = r * dx
      return self._coors

  def to_tuple(self):
      return (self.r1, self.dx, self.jrmt, self.rmt, self.jrws, self.rws)

  def copy(self):
      copy.copy(self)
