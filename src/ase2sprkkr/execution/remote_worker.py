"""Detached worker used by persistent local and remote execution."""

import sys

from .detached import run_detached_worker


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: python -m ase2sprkkr.execution.remote_worker DIRECTORY RUN_ID"
        )
    run_detached_worker(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
