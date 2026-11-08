"""Exceptions for deferred CLI argument validation."""


class CLIUsageError(ValueError):
    """Invalid CLI argument found after lightweight argparse processing."""
