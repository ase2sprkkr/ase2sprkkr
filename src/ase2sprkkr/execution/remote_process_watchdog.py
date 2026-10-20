"""Guard a calculation process group against unexpected supervisor death.

The watchdog leads the calculation's process group within the supervisor's
session. Its only input is a pipe whose writer belongs to the supervisor.
Closing the pipe terminates the group, even when the supervisor gets SIGKILL.
"""

import json
import os
import signal
import sys
import time


def main():
    if sys.argv[1] == "--exec":
        _, _, group_id, error_fd, *argv = sys.argv
        fd = int(error_fd)
        try:
            os.setpgid(0, int(group_id))
            os.set_inheritable(fd, False)
            os.execvp(argv[0], argv)
        except OSError as exc:
            os.write(fd, json.dumps([exc.errno, str(exc)]).encode("utf8"))
            os._exit(127)

    read_fd = int(sys.argv[1])
    os.setpgrp()
    # Survive the first group-wide signal in order to enforce escalation.
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    os.write(sys.stdout.fileno(), b"R")
    try:
        while os.read(read_fd, 65536):
            pass
    finally:
        os.close(read_fd)

    # This watchdog is itself a member and the leader of the target group,
    # so the group cannot disappear while it is sending these signals.
    group = os.getpgrp()
    os.killpg(group, signal.SIGTERM)
    time.sleep(1.0)
    os.killpg(group, signal.SIGKILL)


if __name__ == "__main__":
    main()
