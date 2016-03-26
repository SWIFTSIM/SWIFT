#!/bin/bash
rm brute_force_27_standard.dat swift_dopair_27_standard.dat

./test27cells -p 6 -r 1 -d 0 -f standard

python difffloat.py brute_force_27_standard.dat swift_dopair_27_standard.dat tolerance.dat

exit $?
