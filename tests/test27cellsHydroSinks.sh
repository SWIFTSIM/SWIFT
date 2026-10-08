#!/bin/bash

# Test for the gas-gas neighbour loop used in sink particle formation.
#
# Runs the optimised pair and self interaction functions against the brute-force
# reference and compares the neighbour count (density.wcount) using difffloat.py.
# The comparison must be exact (zero tolerance) since wcount is an integer
# count determined by identical geometric operations on identical input data.

echo "==== test27cellsHydroSinks ===="

# Remove stale output files.
rm -f brute_force_hydro_sinks_27_standard.dat \
       swift_hydro_sinks_dopair_27_standard.dat

echo "Running ./test27cellsHydroSinks -n 5 -r 1 -d 0 -f standard"
./test27cellsHydroSinks -n 5 -r 1 -d 0 -f standard

if [ -e brute_force_hydro_sinks_27_standard.dat ]
then
  if python3 ./difffloat.py brute_force_hydro_sinks_27_standard.dat \
                            swift_hydro_sinks_dopair_27_standard.dat \
                            ./hydro_sinks_tolerance_27_normal.dat 5
  then
    echo "Accuracy test passed"
  else
    echo "Accuracy test FAILED"
    exit 1
  fi
else
  echo "Error: missing output file brute_force_hydro_sinks_27_standard.dat"
  exit 1
fi

echo "------------"

# Perturbed particle positions.
rm -f brute_force_hydro_sinks_27_perturbed.dat \
       swift_hydro_sinks_dopair_27_perturbed.dat

echo "Running ./test27cellsHydroSinks -n 5 -r 1 -d 0.1 -f perturbed"
./test27cellsHydroSinks -n 5 -r 1 -d 0.1 -f perturbed

if [ -e brute_force_hydro_sinks_27_perturbed.dat ]
then
  if python3 ./difffloat.py brute_force_hydro_sinks_27_perturbed.dat \
                            swift_hydro_sinks_dopair_27_perturbed.dat \
                            ./hydro_sinks_tolerance_27_normal.dat 5
  then
    echo "Accuracy test passed"
  else
    echo "Accuracy test FAILED"
    exit 1
  fi
else
  echo "Error: missing output file brute_force_hydro_sinks_27_perturbed.dat"
  exit 1
fi

echo "------------"

exit 0
