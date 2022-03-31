#!/bin/bash
python3 -m pip install `python <<EOF
import configparser
cp=configparser.ConfigParser()
cp.read_file(open('setup.cfg', 'r'))
print(cp['options']['install_requires'])
EOF`
python3 -m pip install pytest
