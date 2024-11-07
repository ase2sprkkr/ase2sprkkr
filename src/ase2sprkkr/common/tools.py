""" This module contain support for the tools """

from pathlib import Path
import pyparsing as pp
import argparse

_unit_registry = None


def parse_inches(string):
    """
    .. doctest::
      >>> parse_inches(1)
      1.0
      >>> int( parse_inches('2cm') * 10000)
      7874
    """
    global _unit_registry
    if isinstance(string, (int, float)):
        return float(string)
    if _unit_registry is None:
        try:
            from pint import UnitRegistry
            _unit_registry = UnitRegistry()
            raise TypeError
        except TypeError:
            result = (pp.Word(pp.nums) + pp.Word(pp.alphas)).parseString(string, parseAll=True)
            out = float(result[0])
            if result[1] == 'inch':
                return out
            if result[1] == 'cm':
                return out / 2.54
            elif result[1]=='mm':
                return out /25.4
            else:
                raise ValueError(f'Pint is not available nad the units {result[1]} is not known')
    out = _unit_registry.parse_expression(string)
    if isinstance(out, (int, float)):
        return float(out)
    return float(out.to(_unit_registry.inch).magnitude)


def parse_tuple_function(type, length=None, max_length=True, delimiter=','):
    """ Returns a function, that can parse a comma delimited tuple of values

    .. doctest ::

    >>> parse_tuple_function(float, 2)("5,4.7")
    (5.0, 4.7)

    >>> parse_tuple_function(float, 3)("1,2")   # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: The given value "1,2" should contain at least 3 values, delimited by ","'

    >>> parse_tuple_function(float, 1)("1,2")   # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: The given value "1,2" should contain no more than 1 values, delimited by ","'

    >>> parse_tuple_function(float, 1,3)("1,2")
    (1.0, 2.0)

    >>> parse_tuple_function(float, 1,1)("1,2")   # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: The given value "1,2" should contain no more than 1 values, delimited by ","'
    """
    if max_length is True:
       max_length = length

    def parse(string):
        out=string.split(delimiter)
        if length and len(out) < length:
           raise ValueError(f'The given value "{string}" should contain at least {length} values, delimited by "{delimiter}"')
        if max_length and len(out) > max_length:
           raise ValueError(f'The given value "{string}" should contain no more than {max_length} values, delimited by "{delimiter}"')
        out = tuple([type(i) for i in out])
        return out

    return parse


def append_id_to_filename(filename, id, connector='_'):
    p = Path(filename)
    return f"{Path.joinpath(p.parent, p.stem)}{connector}{id}{p.suffix}"


string = pp.Regex('"([^"]|"")*"').set_parse_action(lambda x: x[0][1:-1].replace('""','"')) |\
         pp.Regex("'([^']|'')*'").set_parse_action(lambda x: x[0][1:-1].replace("''","'"))
boolean = pp.Keyword('True').set_parse_action(lambda x: True) |\
          pp.Keyword('False').set_parse_action(lambda x: False)


forward = pp.Forward()
token = pp.pyparsing_common.number | pp.Word(pp.alphanums + '-_@#$!/[]') | forward | string | boolean
tupl = pp.Literal('(').suppress() + pp.delimited_list( token, delim = ',' ).set_parse_action(lambda x: tuple(x)) + pp.Literal(')').suppress()
dict_token = ((pp.pyparsing_common.number | pp.Word(pp.alphanums + '-_@#$!/[]')) + pp.Literal(':').suppress() + token).set_parse_action(lambda x: (x[0],x[1]) )
dicti = pp.Literal('{').suppress() + pp.delimited_list( dict_token, delim = ':').set_parse_action(lambda x: dict(x.as_list())) + pp.Literal('}').suppress()
forward << dicti | tupl

option = token ^ pp.Regex('.*').set_parse_action(lambda x: x[0])

equal_value = pp.Literal('=').suppress() + option


def parse_named_option(x:str, numbers_allowed:bool=False):
    """ Parse a given string of the format `name=value`
    If it recognizes number, bool or tuple of values or quoted string in the value,
    it converts it to a given type.

    Returns
    -------
    name:str
      Name of the parsed option
    value:Any
      Value of the parsed option
    """
    if numbers_allowed:
        p = pp.Word(pp.alphanums)
    else:
        p = pp.Word(pp.alphas)
    return tuple((p + equal_value).parse_string(x, True))


def main(local):
    """
    Cli subcommands can be runned on its own.
    This method creates the main function for the the sub-scripts.
    """
    parser = argparse.ArgumentParser(
      description=local['description'],
      formatter_class=argparse.RawDescriptionHelpFormatter
    )
    local['parser'](parser)
    args = parser.parse_args()
    local['run'](args)
