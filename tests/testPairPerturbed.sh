#!/bin/bash
rm brute_force_perturbed.dat swift_dopair_perturbed.dat

./testPairPerturbed -p 6 -r 1

python difffloat.py brute_force_perturbed.dat swift_dopair_perturbed.dat 1e-5 2e-6

exit $?
