""" Formating routines """
import re


def format_for_string(format):
    """
    Return a format string derived from given format,
    that suits for string.

    >>> format_for_string('>14')
    '>14'
    >>> format_for_string('>14e')
    '>14'
    >>> format_for_string('>14e')
    '>14'
    """
    return re.sub('[eEfFgG]$','', format)


def full_format_for_string(format):
    """
    Return a format string derived from the given format,
    that suits for string.

    >>> full_format_for_string('Cokoli {>14} tu')
    'Cokoli {>14} tu'
    >>> full_format_for_string('Remove{>14e}here')
    'Remove{>14}here'
    >>> full_format_for_string('{>14g}')
    '{>14}'
    """
    return re.sub('[eEfFgG]}(?!})','}', format)


def fortran_format(value, format=':.12e'):
    """
    Format the value with a leading zero.

    Format the given float number using the given fortran scientific notation format string,
    but with a leading zero

     >>> fortran_format(0.1)
     '0.1000000000000e-00'

     >>> fortran_format(1, ":> 14.6E")
     '  0.1000000E+01'

    Parameters
    ----------
    value: mixed
      A value to be printed. The value have to be convertable to a float

    format: str
      Fortran format string with a type 'e' or 'E'. Optional

    Return
    ------
    output: str
      A string containing the number in scientific notation with a leading zero.
    """
    la = ('{' + format + '}').format(float(value))
    a = la.lstrip()
    leading = la[:len(la) - len(a)]
    la = a
    a = la.rstrip()
    trailing = la[len(a):]
    if 'E' in format:
        e = a.find('E')
    elif 'e' in format:
        e = a.find('e')
    else:
        raise ("No E in fortran format string: " + format)
    return leading + '0.{}{}{}{:02d}'.format(a[0],a[2:e],a[e:e + 2],abs(int(a[e + 1:]) * 1 + 1)) + trailing
