#!/bin/bash 

#system args:
#1)input_file
#2)output_file
#3)mass degradation slope (linear, quadratic or cubic)
#4)halo mass (H13, H13.5, H14, H15)
#5)halo concentration (see Nobels+2022 and Husko+2023)
#6)central Temp exponent
#7) high resolution region radius
 
python3 make_IC.py M5_lowres8_Tmin7.hdf5 M5_quadratic500_lowres8_Tmin7.hdf5 linear H14 5.6 7 500
