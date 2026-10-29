from tempfile import NamedTemporaryFile
from unittest.mock import patch

from ase.build import bulk
import pytest

from .. import _preparation as preparation
from ..calculator import SPRKKR


class _StopPotentialOpen(RuntimeError):
    pass


def test_save_input_preserves_explicit_false_file_sentinels():
    calculator = SPRKKR(
        atoms=bulk("Fe"),
        input_file=True,
        potential_file=True,
        empty_spheres=False,
    )
    opened_filenames = []

    with NamedTemporaryFile(mode="w+") as input_file:

        def open_file(filename, *args, **kwargs):
            opened_filenames.append(filename)
            if len(opened_filenames) == 1:
                return input_file
            raise _StopPotentialOpen

        with patch.object(preparation, "open_file", side_effect=open_file):
            with pytest.raises(_StopPotentialOpen):
                calculator.save_input(
                    input_file=False,
                    potential_file=False,
                    output_file=False,
                    return_files=True,
                    empty_spheres=False,
                )

    assert opened_filenames == [False, False]
