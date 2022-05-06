################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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
################################################################################

import h5py
import numpy as np
from optparse import OptionParser
from astropy import units
from astropy import constants

def parse_options():

    usage = "usage: %prog [options] file"
    parser = OptionParser(usage=usage)
        
    
    parser.add_option("--lJ",
                      action="store", 
                      dest="lJ",
                      type="float",
                      default = 0.250,
                      help="Jeans wavelength in box size unit")     

    parser.add_option("--rho",
                      action="store", 
                      dest="rho",
                      type="float",
                      default = 0.1,
                      help="Mean gas density in atom/cm3")  

    parser.add_option("--mass",
                      action="store", 
                      dest="mass",
                      type="float",
                      default = 50,
                      help="Gas particle mass in solar mass") 

    parser.add_option("--level",
                      action="store", 
                      dest="level",
                      type="int",
                      default = 5,
                      help="Resolution level: N = (2**l)**3")  


    parser.add_option("-o",
                      action="store", 
                      dest="outputfilename",
                      type="string",
                      default = "box.dat",
                      help="output filename") 

                                            
    (options, args) = parser.parse_args()

    files = args

    return files, options    
    
    
########################################
# main
########################################

files, opt = parse_options()

# define standard units
UnitMass_in_cgs    = 1.989e43    # 10^10 M_sun in grams
UnitLength_in_cgs  = 3.085678e21 # kpc in centimeters
UnitVelocity_in_cgs= 1e5         # km/s in centimeters per second
UnitCurrent_in_cgs = 1           # Amperes
UnitTemp_in_cgs    = 1           # Kelvin
UnitTime_in_cgs    = UnitLength_in_cgs/UnitVelocity_in_cgs

UnitMass           = UnitMass_in_cgs   *units.g
UnitLength         = UnitLength_in_cgs *units.cm
UnitTime           = UnitTime_in_cgs   *units.s
UnitVelocity       = UnitVelocity_in_cgs   *units.cm/units.s



# Number of particles
N = (2**opt.level)**3     # number of particles

# Mean density
rho = opt.rho  # atom/cc
rho = rho*constants.m_p/units.cm**3  

# Gas particle mass
m = opt.mass # in solar mass
m = m*units.Msun

# Gas mass in the box
M = N*m

# Size of the box 
L = (M/rho)**(1/3.)

# Jeans wavelength in box size unit
lJ = opt.lJ    
lJ = lJ*L

# Gravitational constant
G = constants.G

# Jeans wave number
kJ = 2*np.pi/lJ
# Velocity dispersion 
sigma = np.sqrt(4*np.pi*G*rho)/ kJ



print("Number of particles                   : {}".format(N))
print("Equivalent velocity dispertion        : {}".format(sigma.to(units.m/units.s)))

# Convert to code units
m     = m.to(UnitMass).value
L     = L.to(UnitLength).value
rho   = rho.to(UnitMass/UnitLength**3).value
sigma = sigma.to(UnitVelocity).value

# Generate the particles
pos = np.random.random([N, 3])* np.array([L, L, L])
vel  = np.zeros(N)
mass = np.ones(N)*m 
u    = np.ones(N)*sigma**2
ids  = np.arange(N)  
h    = np.ones(N)* 3*L/N**(1/3.)
rho  = np.ones(N)*rho

print("Inter-particle distance (code unit)   : {}".format(L/N**(1/3.)))


# File
fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3


# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"]      = UnitLength_in_cgs
grp.attrs["Unit mass in cgs (U_M)"]        = UnitMass_in_cgs
grp.attrs["Unit time in cgs (U_t)"]        = UnitTime_in_cgs
grp.attrs["Unit current in cgs (U_I)"]     = UnitCurrent_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs


# Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=vel, dtype="f")
grp.create_dataset("Masses", data=mass, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("Densities", data=rho, dtype="f")

fileOutput.close()

