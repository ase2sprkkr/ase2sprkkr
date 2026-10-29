"""Output printing listeners independent of task-specific parsers."""

import re
import sys
from threading import Lock


class OutputListener:
    """Callable receiving (text, stream kind) from a running calculation."""

    def __call__(self, text, kind):
        raise NotImplementedError

    def flush(self):
        """Finish processing any pending output."""


class PrintAll(OutputListener):
    """Print all received output chunks."""

    def __call__(self, text, kind):
        stream = sys.stderr if kind == "err" else sys.stdout
        stream.write(text)
        stream.flush()


class PrintMatching(OutputListener):
    """Print complete lines matching one of the provided regexes."""

    def __init__(self, patterns):
        self.patterns = tuple(patterns)
        self.pending = {"out": "", "err": ""}
        self.lock = Lock()

    def __call__(self, text, kind):
        if kind not in self.pending or not text:
            return
        with self.lock:
            *lines, remainder = (self.pending[kind] + text).split("\n")
            self.pending[kind] = remainder
            stream = sys.stderr if kind == "err" else sys.stdout
            for line in lines:
                line += "\n"
                if any(p.search(line) for p in self.patterns):
                    stream.write(line)
                    stream.flush()

    def flush(self):
        with self.lock:
            for kind, line in self.pending.items():
                if line and any(p.search(line) for p in self.patterns):
                    stream = sys.stderr if kind == "err" else sys.stdout
                    stream.write(line)
                    stream.flush()
                self.pending[kind] = ""


def output_listeners(value):
    """Return (tuple of OutputListener instances, print_info: bool).

    Strings except 'info' are regular expressions; tuples combine modes.
    """
    values = value if isinstance(value, (tuple, list)) else (value,)
    info = all_output = False
    patterns, listeners = [], []
    for item in values:
        if item is False or item is None:
            continue
        if item is True:
            all_output = True
        elif isinstance(item, str):
            if item == "info":
                info = True
            else:
                patterns.append(re.compile(item))
        elif isinstance(item, re.Pattern):
            if isinstance(item.pattern, bytes):
                raise TypeError("print_output requires text regexes")
            patterns.append(item)
        elif isinstance(item, OutputListener):
            listeners.append(item)
        else:
            raise TypeError("Invalid print_output option: expected bool, 'info', regex, or tuple")
    if all_output:
        listeners.append(PrintAll())
    elif patterns:
        listeners.append(PrintMatching(patterns))
    return tuple(listeners), info
