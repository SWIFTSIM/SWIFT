#!/usr/bin/env python3


import numpy as np
from pNbody import Nbody
from astropy import units as u
import argparse



####################################################################
# option parser
####################################################################

description="""Generate a Plummer sphere"""
epilog     ="""
Examples:
--------
"""

parser = argparse.ArgumentParser(description=description,epilog=epilog,formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument("-n",
                  action="store",
                  dest="n",
                  type=int,
                  default=100000,
                  help="number of particles",
                  metavar=" INT")   

parser.add_argument("--regular-grid",
                    action="store_true",
                    default=False,
                    help="Put particles on a regular grid") 

parser.add_argument("-o",
                    action="store",
                    type=str,
                    dest="outputfilename",
                    default='snap.hdf5',
                    help="Name of the output file")  


opt = parser.parse_args()

 
# define units
u_Length   = 1* u.kpc
u_Mass     = 10**10 * u.M_sun
u_Velocity = 1* u.km/u.s
u_Time     = u_Length/u_Velocity
toMsol     = u_Mass.to(u.M_sun).value

outputfile = opt.outputfilename


# 1D : density profile

xmax  = 10000
x     = np.linspace(0,xmax,1000) # kpc
x0    = 5000                     # kpc
rho0  = 1e-8/3                   # 1e10Msol/kpc^3
c     = 1.0
sigma = 500                      # kpc

def fct_Rho(x):
  return rho0 * (c+1-np.tanh((x-x0)/sigma))

def fct_CumMass(x):
  return rho0*( (c+1)*x - sigma*np.log( np.cosh((x-x0)/sigma) / np.cosh(x0/sigma) ))

def fct_normedCumMass(x):
  xmax = max(x)
  return fct_CumMass(x)/fct_CumMass(xmax)

rho  = fct_Rho(x)
Mtot = fct_CumMass(xmax)
m    = fct_normedCumMass(x)


if opt.regular_grid:

  # generate the positions
  n = opt.n
  X = np.linspace(0,1,n)
  xx = X*xmax
  
  # remove the last point (at the border of the domain)
  xx = xx[:-1]
  n = len(xx)

  rho = fct_Rho(xx)
  mass = Mtot*rho/sum(rho)

else:  

  # generate the positions
  n = opt.n
  X = np.linspace(0,1,n)
  xx = np.interp(X,m,x)
   
  # remove the last point (at the border of the domain)
  xx = xx[:-1]
  n = len(xx)

  mass = Mtot*np.ones(n)/n

  
# create the pNbody object
pos  = np.transpose(np.array([xx, np.zeros(n), np.zeros(n)]))
vel  = np.transpose(np.array([np.zeros(n),  np.zeros(n), np.zeros(n)]))

nb = Nbody(status='new',pos=pos,vel=vel,mass=mass,ftype='swift')
# set particles to gas
nb.set_tpe(0)
nb.u_init   = np.ones(nb.nbody)*0.05  
nb.rsp_init = nb.get_rsp_approximation()

  
# add units
nb.UnitLength_in_cm         = u_Length.to(u.cm).value
nb.UnitMass_in_g            = u_Mass.to(u.g).value
nb.UnitVelocity_in_cm_per_s = u_Velocity.to(u.cm/u.s).value
nb.Unit_time_in_cgs         = u_Time.to(u.s).value

nb.hubblefactorcorrection      = False
nb.comovingtoproperconversion  = False
nb.atime                       = 1

nb.boxsize = xmax

nb.rename(outputfile)
nb.write()

# add dimention in the header
import h5py
with h5py.File(outputfile, "a") as f:  
    header = f["Header"]
    header.attrs["Dimension"] = 1  







