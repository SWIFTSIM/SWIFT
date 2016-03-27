#!/bin/bash
rm brute_force_perturbed.dat swift_dopair_perturbed.dat

./testPair -p 6 -r 1 -d 0.1 -f perturbed

python difffloat.py brute_force_perturbed.dat swift_dopair_perturbed.dat tolerance.dat

exit $?
