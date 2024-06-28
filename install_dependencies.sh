#!/bin/bash
python3 -m pip install toml pytest
python3 -m pip install `python list_dependencies.py`
