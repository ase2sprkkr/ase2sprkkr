"""
This module implement simple class that encapsulates the
warning about non-critical errors (or just suspicious data)
in configuration/data files.
"""

import warnings


class DataValidityWarning(UserWarning):
    """
    This Warning should be issued, if there are some invalid data
    or format problems, that should not yield a "hard error" which
    would prevent parsing the data.
    """

    @classmethod
    def warn(cls, out):
        """ Yield the warning """
        warnings.warn(cls(out))
