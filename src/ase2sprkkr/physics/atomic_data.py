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

@cache
def mendeleev_element(atomic_number):
    """Return optional Mendeleev metadata for a real chemical element.

    This is the shared provider for both the Python API and the GUI. The
    heavyweight third-party database is never imported by normal calculations.
    """
    if not 1 <= int(atomic_number) <= 118:
        raise ValueError(f"Mendeleev has no element with atomic number {atomic_number!r}")
    try:
        from mendeleev import element
    except ImportError as exc:
        raise ImportError(
            "Detailed element information requires the optional dependency. "
            'Install with: pip install "ase2sprkkr[mendeleev]"'
        ) from exc
    return element(int(atomic_number))
