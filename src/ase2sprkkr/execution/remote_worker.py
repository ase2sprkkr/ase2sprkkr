"""Detached worker used by :mod:`ase2sprkkr.execution.remote`."""

import sys

from .remote import run_remote_worker


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: python -m ase2sprkkr.execution.remote_worker DIRECTORY")
    run_remote_worker(sys.argv[1])


if __name__ == "__main__":
    main()
