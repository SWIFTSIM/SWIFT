#!/bin/bash

hg clone https://bitbucket.org/rthompson/pygadgetreader

python3 -m venv pygadenv
source pygadenv/bin/activate
pip3 install swiftsimio

cd pygadgetreader

2to3 -nw *.py
2to3 -nw */*.py

python3 setup.py install

cd ..

deactivate
