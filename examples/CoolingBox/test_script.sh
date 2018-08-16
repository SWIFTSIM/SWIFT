#!/bin/bash

max_number() {
    printf "%s\n" "$@" | sort -g | tail -n1
}

# loop over redshift
for z in $(seq 0.01 5.0 25.01); do
  # loop over solar and zero metal abundances
  for solar in {0..1}; do
    # loop over log_10 hydrogen number density
    for nh_exp in $(seq -3 2.0 -1); do
      #change parameters in yml file for calculating explicit ode solution of cooling
      cd ../CoolingRates
      if [ $solar == 1 ]
      then
        sed -i "/InitMetallicity: / s/: \+[[:alnum:],\.,-]\+/: 0.014/g" testCooling.yml
        sed -i "/InitAbundance_Hydrogen: / s/: \+[[:alnum:],\.,-]\+/: 0.70649785/g" testCooling.yml
        sed -i "/InitAbundance_Helium: / s/: \+[[:alnum:],\.,-]\+/: 0.28055534/g" testCooling.yml
        sed -i "/InitAbundance_Carbon: / s/: \+[[:alnum:],\.,-]\+/: 2.0665436e-3/g" testCooling.yml
        sed -i "/InitAbundance_Nitrogen: / s/: \+[[:alnum:],\.,-]\+/: 8.3562563e-4/g" testCooling.yml
        sed -i "/InitAbundance_Oxygen: / s/: \+[[:alnum:],\.,-]\+/: 5.4926244e-3/g" testCooling.yml
        sed -i "/InitAbundance_Neon: / s/: \+[[:alnum:],\.,-]\+/: 1.4144605e-3/g" testCooling.yml
        sed -i "/InitAbundance_Magnesium: / s/: \+[[:alnum:],\.,-]\+/: 5.907064e-4/g" testCooling.yml
        sed -i "/InitAbundance_Silicon: / s/: \+[[:alnum:],\.,-]\+/: 6.825874e-4/g" testCooling.yml
        sed -i "/InitAbundance_Iron: / s/: \+[[:alnum:],\.,-]\+/: 1.1032152e-3/g" testCooling.yml
        sed -i "/CalciumOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.0941736/g" testCooling.yml
        sed -i "/SulphurOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.6054160/g" testCooling.yml
      elif [ $solar == 0 ]
      then
        sed -i "/InitMetallicity: / s/: \+[[:alnum:],\.,-]\+/: 0.0/g" testCooling.yml
        sed -i "/InitAbundance_Hydrogen: / s/: \+[[:alnum:],\.,-]\+/: 0.752/g" testCooling.yml
        sed -i "/InitAbundance_Helium: / s/: \+[[:alnum:],\.,-]\+/: 0.248/g" testCooling.yml
        sed -i "/InitAbundance_Carbon: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/InitAbundance_Nitrogen: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/InitAbundance_Oxygen: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/InitAbundance_Neon: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/InitAbundance_Magnesium: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/InitAbundance_Silicon: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/InitAbundance_Iron: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" testCooling.yml
        sed -i "/CalciumOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.0941736/g" testCooling.yml
        sed -i "/SulphurOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.6054160/g" testCooling.yml
      fi
      rm cooling_*.dat
      ./testCooling -z $z -d $nh_exp
      cd ../CoolingBox
      cp ../CoolingRates/cooling_output.dat ./
  
      #change parameters in coolingBox.yml
      if [ $solar == 1 ]
      then
        sed -i "/InitMetallicity: / s/: \+[[:alnum:],\.,-]\+/: 0.014/g" coolingBox.yml
        sed -i "/InitAbundance_Hydrogen: / s/: \+[[:alnum:],\.,-]\+/: 0.70649785/g" coolingBox.yml
        sed -i "/InitAbundance_Helium: / s/: \+[[:alnum:],\.,-]\+/: 0.28055534/g" coolingBox.yml
        sed -i "/InitAbundance_Carbon: / s/: \+[[:alnum:],\.,-]\+/: 2.0665436e-3/g" coolingBox.yml
        sed -i "/InitAbundance_Nitrogen: / s/: \+[[:alnum:],\.,-]\+/: 8.3562563e-4/g" coolingBox.yml
        sed -i "/InitAbundance_Oxygen: / s/: \+[[:alnum:],\.,-]\+/: 5.4926244e-3/g" coolingBox.yml
        sed -i "/InitAbundance_Neon: / s/: \+[[:alnum:],\.,-]\+/: 1.4144605e-3/g" coolingBox.yml
        sed -i "/InitAbundance_Magnesium: / s/: \+[[:alnum:],\.,-]\+/: 5.907064e-4/g" coolingBox.yml
        sed -i "/InitAbundance_Silicon: / s/: \+[[:alnum:],\.,-]\+/: 6.825874e-4/g" coolingBox.yml
        sed -i "/InitAbundance_Iron: / s/: \+[[:alnum:],\.,-]\+/: 1.1032152e-3/g" coolingBox.yml
        sed -i "/CalciumOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.0941736/g" coolingBox.yml
        sed -i "/SulphurOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.6054160/g" coolingBox.yml
      elif [ $solar == 0 ]
      then
        sed -i "/InitMetallicity: / s/: \+[[:alnum:],\.,-]\+/: 0.0/g" coolingBox.yml
        sed -i "/InitAbundance_Hydrogen: / s/: \+[[:alnum:],\.,-]\+/: 0.752/g" coolingBox.yml
        sed -i "/InitAbundance_Helium: / s/: \+[[:alnum:],\.,-]\+/: 0.248/g" coolingBox.yml
        sed -i "/InitAbundance_Carbon: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/InitAbundance_Nitrogen: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/InitAbundance_Oxygen: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/InitAbundance_Neon: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/InitAbundance_Magnesium: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/InitAbundance_Silicon: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/InitAbundance_Iron: / s/: \+[[:alnum:],\.,-]\+/: 0.000/g" coolingBox.yml
        sed -i "/CalciumOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.0941736/g" coolingBox.yml
        sed -i "/SulphurOverSilicon: / s/: \+[[:alnum:],\.,-]\+/: 0.6054160/g" coolingBox.yml
      fi
        
      # set starting, ending redshift, how often to write to file
      a_begin=$(python -c "print 1.0/(1.0+$z)")
      a_end=$(python -c "print min(1.0, $a_begin*1.0001)")
      first_ouput_a=$a_begin
      delta_a=$(python -c "print 1.0 + ($a_end/$a_begin - 1.0)/500.")
      sed -i "/a_begin: / s/: \+[[:alnum:],\.,-]\+/: $a_begin/g" coolingBox.yml
      sed -i "/a_end: / s/: \+[[:alnum:],\.,-]\+/: $a_end/g" coolingBox.yml
      sed -i "/scale_factor_first: / s/: \+[[:alnum:],\.,-]\+/: $first_ouput_a/g" coolingBox.yml
      sed -i "/delta_time: / s/: \+[[:alnum:],\.,-]\+/: $delta_a/g" coolingBox.yml
  
      # change hydrogen number density
      nh=$(python -c "print 3.555*10.0**($nh_exp+7)/(1.0+$z)**3")
      sed -i "/^rho =/ s/= \S*/= $nh/g" makeIC.py
      for pressure_index in 7; do
        pressure=$(python -c "print 4.5*10.0**($pressure_index + $nh_exp + 3)/(1.0+$z)")
        # change pressure (and hence energy)
        sed -i "/^P =/ s/= \S*/= $pressure/g" makeIC.py
        python makeIC.py
	rm coolingBox_*hdf5
  
        # run cooling box
        ./run.sh
  
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
