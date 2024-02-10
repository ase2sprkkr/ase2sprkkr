class SuspiciousValueWarning(UserWarning):

   def __init__(self, value, message):
       self.value = value
       super().__init__(message)

   def warn(self, stacklevel=1):
       import warnings
       warnings.warn(self, stacklevel=stacklevel + 1)
