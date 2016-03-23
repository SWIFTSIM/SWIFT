#!/bin/bash
rm brute_force.dat swift_dopair.dat

./testPair -p 6 -r 1

python difffloat.py brute_force.dat swift_dopair.dat 1e-5 5e-6

exit $?
