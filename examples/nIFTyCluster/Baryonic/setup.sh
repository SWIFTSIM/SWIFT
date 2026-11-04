#!/bin/bash

# Install swiftsimio and pygadgetreader, needed by convert_to_swift.py, in a
# virtual environment. The original pygadgetreader repository on Bitbucket no
# longer exists; this fork is maintained and works with Python 3.
python3 -m venv pygadenv
source pygadenv/bin/activate
pip3 install swiftsimio git+https://github.com/dnarayanan/pygadgetreader
deactivate
