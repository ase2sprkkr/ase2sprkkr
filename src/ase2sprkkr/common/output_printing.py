"""Selection and live printing of SPR-KKR process output."""

import re
import sys
from threading import Lock


def normalize_print_output(value):
    """Return (show_all, show_scf_info, regexes) for a print_output option.

    Strings other than 'info' are regexes. Tuples (or lists) combine modes.
    Compiled text regexes are accepted by the Python API.
    """
    values = value if isinstance(value, (tuple, list)) else (value,)
    show_all = show_info = False
    patterns = []
    for item in values:
        if item is True:
            show_all = True
        elif item is False or item is None:
            continue
        elif isinstance(item, str):
            if item == "info":
                show_info = True
            else:
                patterns.append(re.compile(item))
        elif isinstance(item, re.Pattern):
            if isinstance(item.pattern, bytes):
                raise TypeError("print_output requires a text regex, not a bytes regex")
            patterns.append(item)
        else:
            raise TypeError(
                "print_output must be bool, 'info', a regex, or a tuple of these"
            )
    return show_all, show_info, tuple(patterns)


class OutputPrinter:
    """Observer that mirrors complete lines matching print_output regexes.

    Execution delivers arbitrary text chunks, not necessarily whole lines.
    Each output stream is buffered separately and the final unterminated line
    is emitted by flush(). True streams every chunk unchanged.
    """

    def __init__(self, mode):
        self.show_all, self.show_info, self.patterns = normalize_print_output(mode)
        self._pending = {"out": "", "err": ""}
        self._lock = Lock()

    @property
    def enabled(self):
        return self.show_all or bool(self.patterns)

    @staticmethod
    def _stream(kind):
        return sys.stderr if kind == "err" else sys.stdout

    def __call__(self, text, kind):
        if kind not in self._pending or not text:
            return
        with self._lock:
            stream = self._stream(kind)
            if self.show_all:
                stream.write(text)
            elif self.patterns:
                lines = (self._pending[kind] + text).split("\n")
                self._pending[kind] = lines.pop()
                for line in lines:
                    line += "\n"
                    if any(pattern.search(line) for pattern in self.patterns):
                        stream.write(line)

    def flush(self):
        """Emit a matching last line without a terminating newline, once."""
        with self._lock:
            if not self.patterns or self.show_all:
                return
            for kind, line in self._pending.items():
                if line and any(pattern.search(line) for pattern in self.patterns):
                    self._stream(kind).write(line)
                self._pending[kind] = ""
