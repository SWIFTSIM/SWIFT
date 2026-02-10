#!/usr/bin/env python3


import numpy as np
from pNbody import ic
from astropy import units as u

# define units
u_Length   = 1* u.kpc
u_Mass     = 10**10 * u.M_sun
u_Velocity = 1* u.km/u.s
u_Time     = u_Length/u_Velocity
toMsol     = u_Mass.to(u.M_sun).value
G          = 4.299581e+04

M = 0.453    # 1e10 Msol
R = 0.3      # kpc
L = 10       # boxsize [kpc}

V   = 4/3.*np.pi*R**3
rho = M/V
tff = np.sqrt(3*np.pi/(32*G*rho))
print(tff)

n = 100000
outputfile = "snap.hdf5"


# create the pNbody object
nb = ic.homosphere(n,R,R,R,irand=1,ftype='swift')
nb.mass = np.ones(n)*M/n

hydro=True
if hydro:
  # set particles to gas
  nb.set_tpe(0)
  nb.u_init   = np.ones(nb.nbody)*0.05  
  nb.rsp_init = nb.get_rsp_approximation()
else:
  # set particles to dm
  nb.set_tpe(1)
  

# add units
nb.UnitLength_in_cm         = u_Length.to(u.cm).value
nb.UnitMass_in_g            = u_Mass.to(u.g).value
nb.UnitVelocity_in_cm_per_s = u_Velocity.to(u.cm/u.s).value
nb.Unit_time_in_cgs         = u_Time.to(u.s).value

nb.hubblefactorcorrection      = False
nb.comovingtoproperconversion  = False
nb.atime                       = 1

nb.boxsize = np.array([L,L,L])

nb.rename(outputfile)
nb.write()





