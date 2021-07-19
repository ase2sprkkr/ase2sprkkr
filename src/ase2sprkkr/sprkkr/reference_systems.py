import copy

class ReferenceSystem:
  """
  #TODO - doc
  """

  def __init__(self, vref, rmtref):
      self.vref = vref
      self.rmtref = rmtref

  def to_tuple(self):
      return self.vref, self.rmtref

  @staticmethod
  def default():
      return ReferenceSystem(4.0,0.0)

  def copy(self):
      return copy.copy(self)
