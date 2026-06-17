#!/usr/bin/env python3

import agama
import numpy as np
from pNbody import ic


# number of particles
n = 10000

# Plummer scale radius
Rc   = 0.1  # in kpc

# Maximal radius
Rmax = 20   # in kpc

# Total Mass
Mtot = 0.01 # in 1e10 Msol

# create the static nbody model
nb = ic.plummer(n,1,1,1,Rc,Rmax,Mtot)

# Define the physical units used in the code: the choice below corresponds to
# length scale = 1 kpc, velocity = 1 km/s, mass = 1e10 Msun
agama.setUnits(mass=1e10,length=1,velocity=1)

# Create a spherical potential based on the particle distribution
potential = agama.Potential(type='multipole',particles=(nb.pos, nb.mass),lmax=4, symmetry='s',rmin=0., rmax=Rmax)

potential.export("potential.txt")
print('Done, enjoy your potential!')
