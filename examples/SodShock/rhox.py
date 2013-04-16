###############################################################################
 # This file is part of SWIFT.
 # Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 # 
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU Lesser General Public License as published
 # by the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 # 
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 # 
 # You should have received a copy of the GNU Lesser General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 # 
 ##############################################################################

import h5py
import random
import sys
import math
from numpy import *

# Reads the HDF5 output of SWIFT and generates a radial density profile
# of the different physical values.

# Input values?
if len(sys.argv) < 3 :
    print "Useage: " , sys.argv[0] , " <filename> <nr. bins>"
    exit()
    
# Get the input arguments
fileName = sys.argv[1];
nr_bins = int( sys.argv[2] );


# Open the file
fileName = sys.argv[1];
file = h5py.File( fileName , 'r' )

# Get the space dimensions.
grp = file[ "/Header" ]
boxsize = grp.attrs[ 'BoxSize' ]
boxsize = boxsize[0]

# Get the particle data
grp = file.get( '/PartType0' )
ds = grp.get( 'Coordinates' )
coords = ds[...]
ds = grp.get( 'Velocities' )
v = ds[...]
# ds = grp.get( 'Mass' )
# m = ds[...]
ds = grp.get( 'SmoothingLength' )
h = ds[...]
ds = grp.get( 'InternalEnergy' )
u = ds[...]
ds = grp.get( 'ParticleIDs' )
ids = ds[...]
ds = grp.get( 'Density' )
rho = ds[...]

# Get the maximum radius
r_max = boxsize

# Init the bins
nr_parts = coords.shape[0]
bins_v = zeros( nr_bins )
bins_m = zeros( nr_bins )
bins_h = zeros( nr_bins )
bins_u = zeros( nr_bins )
bins_rho = zeros( nr_bins )
bins_count = zeros( nr_bins )
bins_P = zeros( nr_bins )

# Loop over the particles and fill the bins.
for i in range( nr_parts ):

    # Get the box index.
    r = coords[i,0]
    ind = floor( r / r_max * nr_bins )
    
    # Update the bins
    bins_count[ind] += 1
    bins_v[ind] += v[i,0] # sqrt( v[i,0]*v[i,0] + v[i,1]*v[i,1] + v[i,2]*v[i,2] )
    # bins_m[ind] += m[i]
    bins_h[ind] += h[i]
    bins_u[ind] += u[i]
    bins_rho[ind] += rho[i]
    bins_P[ind] += (2.0/3)*u[i]*rho[i]
    
# Loop over the bins and dump them
print "# bucket left right count v m h u rho"
for i in range( nr_bins ):

    # Normalize by the bin volume.
    r_left = r_max * i / nr_bins
    r_right = r_max * (i+1) / nr_bins
    vol = 4/3*math.pi*(r_right*r_right*r_right - r_left*r_left*r_left)
    ivol = 1.0 / vol

    print "%i %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e" % \
        ( i , r_left , r_right , \
          bins_count[i] * ivol , \
          bins_v[i] / bins_count[i] , \
          bins_m[i] * ivol , \
          bins_h[i] / bins_count[i] , \
          bins_u[i] / bins_count[i] , \
          bins_rho[i] / bins_count[i] ,
          bins_P[i] / bins_count[i] )
    
    
