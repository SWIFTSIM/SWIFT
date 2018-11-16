#!/bin/bash

# loop over redshift
for z in 0.1; do
  # loop over solar and zero metal abundances
  for solar in 1; do
    # loop over log_10 hydrogen number density
    for nh_exp in -1; do
      # Choose which yml file to use for different metallicities
      if [ $solar == 1 ]
      then
        yml_file="coolingBox_solar.yml"
      elif [ $solar == 0 ]
      then
        yml_file="coolingBox_primordial.yml"
      fi

      # Calculate cooling rates for given redshift, density and metallicity
      rm cooling_*.dat
      ./testCooling -z $z -d $nh_exp -m $yml_file
  
        
      # set starting, ending redshift, how often to write to file
      a_begin=$(python -c "print 1.0/(1.0+$z)")
      a_end=$(python -c "print min(1.0, $a_begin*1.0001)")
      delta_a=$(python -c "print 1.0 + ($a_end/$a_begin - 1.0)/50.")
  
      # change hydrogen number density
      if [ $solar == 0 ]; 
      then
        rho=$(python -c "print 10.0**$nh_exp/0.7*1.6726e-24")
      else
        rho=$(python -c "print 10.0**$nh_exp/0.752*1.6726e-24")
      fi
      for pressure_index in 8; do
        pressure=$(python -c "print 6.68e-14")
        
	python makeIC.py $rho $pressure
	rm coolingBox_*hdf5
  
        # run cooling box
	../swift -s -c -C -t 4 -P Cosmology:a_begin:$a_begin -P Cosmology:a_end:$a_end -P Snapshots:scale_factor_first:$a_begin -P Snapshots:delta_time:$delta_a -P Statistics:scale_factor_first:$a_begin -P Statistics:delta_time:$delta_a $yml_file 2>&1 | tee output.log
  
        max=0
        for file in $( ls -lth coolingBox_*.hdf5 | head -n 1 ); do
          f=${file//hdf5/}
          f=${f//[!0-9]/}
        done
        frames=$((10#$f))
	echo "number of frames $frames"

        # check if everything worked and create plots
        python analytical_test.py $z $nh_exp $pressure_index $solar $frames
      done
    done
  done
done
