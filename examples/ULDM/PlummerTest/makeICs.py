#!/usr/bin/env python3


import numpy as np
from pNbody import ic
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
G          = 4.299581e+04

M = 1        # 1e10 Msol
a = 0.1      # scale radius
R = 50       # kpc
Rt = 5       # truncation radius
L = 20       # boxsize [kpc}

outputfile = opt.outputfilename


# create the pNbody object
nb = ic.plummer(opt.n,1,1,1,a,R,irand=1,ftype='swift')
nb.mass = np.ones(nb.nbody)*M/nb.nbody
nb = nb.selectc(nb.rxyz()<Rt)

hydro=True
if hydro:
  # set particles to gas
  nb.set_tpe(0)
  nb.u_init   = np.ones(nb.nbody)*0.05  
  nb.rsp_init = nb.get_rsp_approximation()
  #nb.rsp_init = np.clip(nb.rsp_init,0,5)
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





