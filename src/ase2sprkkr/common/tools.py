""" This module contain support for the tools """

from pathlib import Path
import pyparsing as pp

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

string = pp.Regex('"([^"]|"")*"').set_parse_action(lambda x: x[0][1:-1].replace('""','"')) |\
         pp.Regex("'([^']|'')*'").set_parse_action(lambda x: x[0][1:-1].replace("''","'"))
boolean = pp.Keyword('True').set_parse_action(lambda x: True) |\
          pp.Keyword('False').set_parse_action(lambda x: False)
token = pp.pyparsing_common.number | pp.Word(pp.alphanums + '-_@#$!/[]') | string | boolean
tupl =  pp.Literal('(').suppress() + pp.delimited_list( token, delim = ',' ).set_parse_action(lambda x: tuple(x)) + pp.Literal(')').suppress()
option = (token | tupl) ^ pp.Regex('.*').set_parse_action(lambda x: x[0])

name_value = pp.Word(pp.alphas) + pp.Literal('=').suppress() + option


def parse_named_option(x:str):
    """ Parse a given string of the format `name=value`
    If it recognize number, bool or tuple of values or quoted string in the value,
    convert it to a given type.

    Returns
    -------
    name:str
      Name of the parsed option
    value:Any
      Value of the parsed option
    """
    return tuple(name_value.parse_string(x, True))
