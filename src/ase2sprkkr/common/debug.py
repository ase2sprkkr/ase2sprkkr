import sys


def add_debug_hook(only_original=False):
    """
    Breaks when an unhandled exception occurs, if
    the terminal is present

    https://stackoverflow.com/questions/1237379/how-do-i-set-sys-excepthook-to-invoke-pdb-globally-in-python
    """

    def info(type, value, tb):
        if (  # hasattr(sys, "ps1") or
            not sys.stderr.isatty() or not sys.stdin.isatty()
        ):
            # stdin or stderr is redirected, just do the normal thing
            original_hook(type, value, tb)
        else:
            # a terminal is attached and stderr is not redirected, debug
            import traceback
            import pdb

            traceback.print_exception(type, value, tb)
            print
            pdb.pm()

    original_hook = sys.excepthook
    if not only_original or sys.excepthook == sys.__excepthook__:
        # if someone already patched excepthook, let them win
        sys.excepthook = info


def backtrace_pyparsing(term, fn):
    fn(term)
    if hasattr(term, "exprs"):
        for i in term.exprs:
            debug_pyparsing(i)
    else:
        term = getattr(term, "expr", None)
        if term:
            debug_pyparsing(term)


def debug_pyparsing(term):
    backtrace_pyparsing(term, lambda term: term.set_debug(True))


def check_whitespaces(term, chars={'\t', ' ', '\r'}):

    def check(term):
        if term.whiteChars != chars:
            raise ValueError('Bad white chars')
        backtrace_pyparsing(term, check)


def print_grammar(expr, indent=0, stack=None):
    if stack is None:
        stack = set()

    prefix = "  " * indent

    if id(expr) in stack:
        print(f"{prefix}{type(expr).__name__}: {expr}  <recursive>")
        return

    print(f"{prefix}{type(expr).__name__}: {expr}")

    stack.add(id(expr))
    try:
        if hasattr(expr, "exprs"):
            for child in expr.exprs:
                print_grammar(child, indent + 1, stack)

        elif hasattr(expr, "expr"):
            if expr.expr is not None:
                print_grammar(expr.expr, indent + 1, stack)
    finally:
        stack.remove(id(expr))
