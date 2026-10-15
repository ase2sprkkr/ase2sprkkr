"""Lightweight access to elemental properties not included in ASE."""

from functools import cache


@cache
def periodic_table():
    """Return the shared Elementy periodic table."""
    from elementy import PeriodicTable

    return PeriodicTable()


def valence_electrons(symbol):
    """Return the conventional number of valence electrons for an element."""
    try:
        return int(periodic_table().elements[symbol]["valence_electrons"])
    except KeyError as exc:
        raise ValueError(f"Unknown chemical element: {symbol!r}") from exc
