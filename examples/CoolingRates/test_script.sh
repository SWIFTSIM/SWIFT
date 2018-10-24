#!/bin/bash

#redshift_array=( 0 0.049 0.101 0.155 0.211 0.271 0.333 0.399 0.468 0.540)
redshift_array=( 1 2 3 4 5 6 7 8 9 10 )

for i in $(seq 10 -1 0); do
  redshift=${redshift_array[$i]}
  ./testCooling -z $redshift -d 3
  file='cooling_output_'$i'.dat'
  cp cooling_output.dat $file
done

#for nh in $(seq -1 -1 -6); do
#  ./testCooling -z 1 -d $nh
#  file='cooling_output_'$(expr 0 - $nh)'.dat'
#  cp cooling_output.dat $file
#done

python cooling_rates_plot.py
