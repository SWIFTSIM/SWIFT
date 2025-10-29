#!/bin/bash

# make current swift executable clean
make clean

# configure with desired options
./configure --disable-mpi \
            --with-hydro=sphenix \
            --with-kernel=wendland-C2 \
            --with-hydro-dimension=3 \
            --with-ext-potential=isothermal 
            

# build            
make