"""Helpers for scoped temporary configuration values."""

from contextlib import contextmanager


@contextmanager
def temporary_option_value(option, value):
    """Temporarily override one option's effective value.

    The override is converted and validated like a normal option value, but is
    stored in the transient ``result`` slot so the canonical stored value is not
    changed. The exact original option state, including an existing transient
    result, is restored even if the body raises.

    Unlike :class:`ConfigurationTransaction`, this helper is deliberately
    option-local. A configuration transaction is root-scoped and would also
    roll back unrelated mutations made to the same configuration while the
    temporary value is active.
    """
    change = option.Change(option)
    option.result = option._pack_value(value)
    try:
        yield option
    finally:
        change.rollback()
