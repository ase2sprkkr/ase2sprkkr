import sys
from pathlib import Path

def patch_package():
    file = Path(__file__).resolve()
    top = file.parents[3]
    sys.path.append(str(top))
    return 'sprkkr.sprkkr_io.tests'

