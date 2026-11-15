"""Output printing listeners independent of task-specific parsers."""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from threading import Lock


class OutputListener:
    """Callable receiving (text, stream kind) from a running calculation."""

    def __call__(self, text, kind):
        raise NotImplementedError

    def flush(self):
        """Finish processing any pending output."""


@dataclass(frozen=True)
class NormalizedOutputFilter:
    """Resolved output modes, independent of how the user specified them.

    A single regex is kept as a compiled Pattern, several regexes as a
    tuple of compiled Patterns. 'info' is always a separate control flag.
    """

    all_output: bool = False
    info: bool = False
    regex: re.Pattern | tuple[re.Pattern, ...] | None = None
    listeners: tuple[OutputListener, ...] = ()

    @property
    def regexes(self) -> tuple[re.Pattern, ...]:
        if self.regex is None:
            return ()
        if isinstance(self.regex, re.Pattern):
            return (self.regex,)
        return self.regex

    @property
    def matcher(self):
        """Return a reusable predicate, using search rather than match."""
        if self.all_output:
            return lambda line: True
        patterns = self.regexes
        if not patterns:
            return lambda line: False
        if len(patterns) == 1:
            search = patterns[0].search
            return lambda line: search(line) is not None
        return lambda line: any(pattern.search(line) for pattern in patterns)

    @property
    def mode(self) -> str:
        """Mode for a simple GUI; custom listeners have no GUI equivalent."""
        if self.all_output:
            return "all"
        if self.regex is not None:
            return "info_regex" if self.info else "regex"
        return "info_regex" if self.info else "none"

    @property
    def regex_text(self) -> str:
        return serialize_output_regexes(self.regexes)


def normalize_output_filter(value) -> NormalizedOutputFilter:
    """Normalize the public print_output argument without merging regex ASTs.

    'info' is a reserved string (never a regex); a compiled regex named
    'info' remains a regex. Explicit True makes output regexes and the
    summary redundant. Custom output listeners retain their own behavior.
    """
    if isinstance(value, NormalizedOutputFilter):
        return value
    items = value if isinstance(value, (tuple, list)) else (value,)
    info = all_output = False
    patterns = []
    listeners = []
    seen = set()
    for item in items:
        if item is False or item is None:
            continue
        if item is True:
            all_output = True
        elif isinstance(item, str):
            if item == "info":
                info = True
                continue
            pattern = re.compile(item)
            key = (pattern.pattern, pattern.flags)
            if key not in seen:
                seen.add(key)
                patterns.append(pattern)
        elif isinstance(item, re.Pattern):
            if isinstance(item.pattern, bytes):
                raise TypeError("print_output requires text regexes")
            key = (item.pattern, item.flags)
            if key not in seen:
                seen.add(key)
                patterns.append(item)
        elif isinstance(item, OutputListener):
            listeners.append(item)
        else:
            raise TypeError(
                "Invalid print_output option: expected bool, 'info', regex, or tuple"
            )
    if all_output:
        # Every raw line is already printed, so no matching regex is needed.
        # The GUI deliberately treats all + summary as simply all.
        return NormalizedOutputFilter(True, False, None, tuple(listeners))
    regex = (patterns[0] if len(patterns) == 1
             else tuple(patterns) if patterns else None)
    return NormalizedOutputFilter(False, info, regex, tuple(listeners))


def _pattern_source(pattern: re.Pattern) -> str:
    """Add flags passed to re.compile() that aren't already in its source."""
    source = pattern.pattern
    inferred = re.compile(source).flags
    extra = pattern.flags & ~inferred
    available = ((re.I, "i"), (re.M, "m"), (re.S, "s"),
                 (re.X, "x"), (re.A, "a"))
    known = 0
    flags = ""
    for bit, letter in available:
        known |= bit
        if extra & bit:
            flags += letter
    if extra & ~(known | re.DEBUG):
        raise ValueError("Regex has flags that cannot be serialized")
    if flags:
        source = f"(?{flags}){source}"
    return source


def serialize_output_regexes(patterns: tuple[re.Pattern, ...]) -> str:
    """Encode regexes for one text field, distinguishing one from many.

    Several patterns are written as "pattern1"|"pattern2"; doubled quotes
    inside quoted patterns represent literal quotes. A single pattern
    beginning with a double quote gets a regex backslash escape.
    """
    sources = tuple(_pattern_source(pattern) for pattern in patterns)
    if not sources:
        return ""
    if len(sources) == 1 and sources[0] and not sources[0].startswith('"'):
        return sources[0]
    if len(sources) == 1 and sources[0].startswith('"'):
        return '\\' + sources[0]
    return "|".join('"' + source.replace('"', '""') + '"'
                    for source in sources)


def parse_output_regexes(text: str) -> tuple[re.Pattern, ...]:
    """Read serialized GUI regexes, or compile a single unquoted regex."""
    if not text.startswith('"'):
        return (re.compile(text),) if text else ()
    patterns = []
    pos = 0
    while pos < len(text):
        if text[pos] != '"':
            raise ValueError("Expected a quoted regular expression")
        pos += 1
        value = []
        while True:
            if pos == len(text):
                raise ValueError("Unterminated quoted regular expression")
            char = text[pos]
            pos += 1
            if char != '"':
                value.append(char)
            elif pos < len(text) and text[pos] == '"':
                value.append('"')
                pos += 1
            else:
                break
        patterns.append(re.compile("".join(value)))
        if pos == len(text):
            break
        if text[pos] != "|":
            raise ValueError("Expected | between quoted regular expressions")
        pos += 1
        if pos == len(text):
            raise ValueError("Missing regular expression after |")
    return tuple(patterns)


def output_filter_from_gui(mode: str, regex_text: str):
    """Convert the simple GUI's mode and regex syntax to a config value."""
    if mode == "all":
        return True
    if mode == "none":
        return False
    if mode not in ("regex", "info_regex"):
        raise ValueError(f"Unknown output mode: {mode}")
    if not regex_text.strip() and mode == "info_regex":
        return "info"
    if not regex_text.strip():
        raise ValueError("Enter a regular expression for the selected output mode")
    try:
        patterns = parse_output_regexes(regex_text)
    except re.error as exc:
        raise ValueError(f"Invalid output regular expression: {exc}") from exc
    # String 'info' denotes a summary, not a regex. Wrap that particular
    # regex in a neutral noncapturing group when serializing to config.
    sources = tuple(
        "(?:info)" if (source := _pattern_source(p)) == "info" else source
        for p in patterns
    )
    if mode == "info_regex":
        return ("info", *sources)
    return sources[0] if len(sources) == 1 else sources


class PrintAll(OutputListener):
    """Print all received output chunks."""

    def __call__(self, text, kind):
        stream = sys.stderr if kind == "err" else sys.stdout
        stream.write(text)
        stream.flush()


class PrintMatching(OutputListener):
    """Print complete lines matching one of the provided regexes."""

    def __init__(self, patterns):
        normalized = normalize_output_filter(tuple(patterns))
        self.patterns = normalized.regexes
        self._matches = normalized.matcher
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
                if self._matches(line):
                    stream.write(line)
                    stream.flush()

    def flush(self):
        with self.lock:
            for kind, line in self.pending.items():
                if line and self._matches(line):
                    stream = sys.stderr if kind == "err" else sys.stdout
                    stream.write(line)
                    stream.flush()
                self.pending[kind] = ""



def output_listeners(value):
    """Return (tuple of OutputListener instances, print_info: bool).

    The same normalization is also used by GUI configuration editors.
    """
    normalized = normalize_output_filter(value)
    listeners = list(normalized.listeners)
    if normalized.all_output:
        listeners.append(PrintAll())
    elif normalized.regex is not None:
        listeners.append(PrintMatching(normalized.regexes))
    return tuple(listeners), normalized.info
