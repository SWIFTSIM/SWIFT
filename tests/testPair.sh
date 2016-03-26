#!/bin/bash
rm brute_force_standard.dat swift_dopair_standard.dat

./testPair -p 6 -r 1 -d 0 -f standard

python difffloat.py brute_force_standard.dat swift_dopair_standard.dat tolerance.dat

exit $?
