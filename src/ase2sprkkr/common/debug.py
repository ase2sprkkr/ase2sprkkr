import sys


def add_debug_hook(only_original=False):
    """
    Breaks when an unhandled exception occurs, if
    the terminal is present

    https://stackoverflow.com/questions/1237379/how-do-i-set-sys-excepthook-to-invoke-pdb-globally-in-python
    """

    def info(type, value, tb):
       if (  # hasattr(sys, "ps1") or
          not sys.stderr.isatty() or
          not sys.stdin.isatty()):
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


def debug_pyparsing(term):
    term.setDebug(True)
    if hasattr(term, 'exprs'):
        for i in term.exprs:
            debug_pyparsing(i)
    elif hasattr(term, 'expr'):
        debug_pyparsing(term.expr)
