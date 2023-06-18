""" This module contain support for the tools """

from pathlib import Path

def parse_tuple_function(type, length=None, max_length=None, delimiter=','):
    """ Returns a function, that can parse a comma delimited tuple of values """
    if max_length is None:
       max_length = length
    def parse(string):
        out=string.split(delimiter)
        if length and len(out) < length:
           raise ValueException(f"The given value '{string}' should contain at least {length} values, delimited by '{delimiter}")
        if max_length and len(out) > length:
           raise ValueException(f"The given value '{string}' should contain no more than {max_length} values, delimited by '{delimiter}")
        out = tuple([type(i) for i in out])
        return out


def append_id_to_filename(filename, id, connector='_'):
    p = Path(filename)
    return f"{Path.joinpath(p.parent, p.stem)}{p.stem}{connector}{id}{p.suffix}"
