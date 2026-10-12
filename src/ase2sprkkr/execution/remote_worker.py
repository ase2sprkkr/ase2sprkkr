"""Detached worker used by :mod:`ase2sprkkr.execution.remote`."""

import sys

from .remote import run_remote_worker


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python -m ase2sprkkr.execution.remote_worker DIRECTORY RUN_ID"
        )
    run_remote_worker(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
