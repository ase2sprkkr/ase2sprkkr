class Mesh:
  @staticmethod
  def default():
      """ SPR-KKR computes the mesh itself, if zeros are given """
      return ExponentialMesh(0.,0.,0,0.,0,0.)

class ExponentialMesh(Mesh):
  """ 
  Radial mesh definition for an atomic site 
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
      self.r1 = r1
      self.dx = dx
      self.jrmt = jrmt
      self.rmt = rmt
      self.jrws = jrws
      self.rws = rws

  def to_tuple(self):
      return (self.r1, self.dx, self.jrmt, self.rmt, self.jrws, self.rws)