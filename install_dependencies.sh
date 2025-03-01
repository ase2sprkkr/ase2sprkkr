#!/bin/bash
python3 -m pip install git+https://github.com/lokik/spglib.git
python3 -m pip install . --verbose
python3 << EOF
import numpy
print(numpy.__version__)
print(numpy.__path__)
EOF
