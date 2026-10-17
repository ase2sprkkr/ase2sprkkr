"""Entry point for the crash-surviving local detached helper."""

from __future__ import annotations

import sys

from .detached import run_detached_worker


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if len(argv) != 2:
        raise SystemExit("usage: detached_worker DIRECTORY RUN_ID")
    run_detached_worker(argv[0], argv[1])


if __name__ == "__main__":
    main()
