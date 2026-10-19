"""Guard a calculation process group against unexpected supervisor death.

The watchdog leads the calculation's process group within the supervisor's
session. Its only input is a pipe whose writer belongs to the supervisor.
Closing the pipe terminates the group, even when the supervisor gets SIGKILL.
"""

import os
import signal
import sys
import time


def main():
    read_fd = int(sys.argv[1])
    # Survive the first group-wide signal in order to enforce escalation.
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    os.write(sys.stdout.fileno(), b"R")
    try:
        while os.read(read_fd, 65536):
            pass
    finally:
        os.close(read_fd)

    group = os.getpgrp()
    try:
        os.killpg(group, signal.SIGTERM)
    except ProcessLookupError:
        return
    time.sleep(1.0)
    try:
        os.killpg(group, signal.SIGKILL)
    except ProcessLookupError:
        pass


if __name__ == "__main__":
    main()
