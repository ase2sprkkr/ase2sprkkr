import sys

import pytest


_DEFAULT_CALC_ARGS = {
    "input_file": "output_test_calc.inp",
    "output_file": "output_test_calc.out",
    "potential_file": "output_test_calc.pot",
    "print_output": "-v" in sys.argv or "--verbose" in sys.argv,
    "mpi": "auto",
    "options": {"NKTAB": 5, "NE": 20},
    "empty_spheres": False,
}


@pytest.fixture
def calc_args(tmp_path):
    """Build calculator arguments in an isolated pytest temporary directory."""

    def make(kind=None, **overrides):
        args = dict(_DEFAULT_CALC_ARGS)
        args["directory"] = str(tmp_path)
        args["options"] = dict(_DEFAULT_CALC_ARGS["options"])

        options = overrides.pop("options", None)
        if options is not None:
            args["options"].update(options)
        args.update(overrides)

        if kind == "save_input":
            args.pop("print_output", None)
            args.pop("mpi", None)
        return args

    return make
