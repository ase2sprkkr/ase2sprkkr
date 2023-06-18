""" Routines for dealing with symbols """
from ase.atoms import Symbols
import re

def filename_from_symbols(symbols, max_len=10):
    """
    Choose a best short representation of given Symbols.

    >>> filename_from_symbols("GeXTeXGeXTeXGeXTeXGeXTeX9", max_len=12)
    '{GeXTeX}4X8'
    >>> filename_from_symbols("GeXTeXGeXTeXGeXTeXGeXTeX9")
    'Ge4Te4X16'
    >>> filename_from_symbols("GeXTeXGeXTeXH")
    '{GeXTeX}2H'
    >>> filename_from_symbols("GeXTeXGeXTeX")
    'GeXTeX'
    >>> filename_from_symbols("C4H4OC4H4OC2C4H4OC4H4OC2")
    'C20H16O4'
    >>> filename_from_symbols("CO2")
    'CO2'
    """
    if not isinstance(symbols, Symbols):
        symbols = Symbols.fromsymbols(symbols)
    out = symbols.get_chemical_formula('reduce')
    if len(out) <= max_len:
       return out
    out = pretty_symbols(symbols)
    out = re.sub(r"^\{([^}]*)\}[0-9]+$", r"\1", out)
    if len(out) <= max_len:
       return out
    out = symbols.get_chemical_formula()
    if len(out) <= max_len:
       return out
    return symbols.get_chemical_formula(empirical=True)

def pretty_symbols(symbols):
    """
    Make a symbols string more pretty. Especially suitable for a long 2D semiinfinite bulks.

    >>> pretty_symbols("GeXTeXGeXTeXGeXTeXGeXTeX9")
    '{GeXTeX}4X8'
    >>> pretty_symbols("CO2")
    'CO2'

    #TODO - this test does not work yet:
    > >> pretty_symbols("C4H4OC4H4OC2C4H4OC4H4OC2")
    '{{C4H4O}2C2}2'
    """
    symbols=re.sub("([A-Z][a-z]*)([0-9]+)", lambda m: m.group(1)*int(m.group(2)), str(symbols))

    prev=[]
    i = 0
    while i < len(symbols):
        for t in prev:
            lt = len(t)
            j = i+lt
            if t == symbols[i:j]:
               k = j+lt
               while t == symbols[j:k]:
                     j=k
                     k+=lt
               repeat = (j - i) // lt + 1
               t = pretty_symbols(t)
               if re.match('^[A-Z][a-z]*$', t):
                   sub = f'{t}{repeat}'
               else:
                   sub = f'{{{t}}}{repeat}'
               symbols = f'{symbols[:i-lt]}{sub}{symbols[j:]}'
               i = i - lt + len(sub)
               prev = [ p[:-lt] + sub for p in prev if len(p) >= lt ]
               break
        else:
               prev = [ p + symbols[i] for p in prev ]
               prev.append(symbols[i])
               i=i+1
    return symbols
