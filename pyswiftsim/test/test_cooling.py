#!/usr/bin/env python3
from pyswiftsim import pointer
from pyswiftsim import cooling

import astropy.units

import matplotlib.pyplot as plt

import numpy as np

# number of particles
N  = 100
# params file
params = "/home/loikki/swift_test/cooling_box/coolingBox.yml"

# read params
params = cooling.read_params(params)
#print(params)

# init units
d = cooling.init_units(params)
units = d["units"]
consts = d["constants"]

# init cooling
cooling_data = cooling.pycooling_init(params, units, consts)

# Init variables
rho = np.logspace(-6, 4, N).astype(np.float) # in cm-3
rho *= consts.const_proton_mass * units.UnitLength_in_cgs**3
T = np.logspace(1,9, N).astype(np.float) # in K

rho, T = np.meshgrid(rho, T)
shape = rho.shape
rho = rho.reshape(np.prod(shape))
T = T.reshape(np.prod(shape))

# u = k_b T / (gamma - 1) mu m_p
mu = 2
u = consts.const_boltzmann_k * T / ((cooling._hydro_gamma - 1.) * mu * consts.const_proton_mass)

# compute cooling
rate = cooling.pycooling_rate(units, cooling_data, consts, rho, u)

cs = cooling.soundspeed_from_internal_energy(rho, u)

rho = rho.reshape(shape) / (units.UnitLength_in_cgs**3 * consts.const_proton_mass)
T = T.reshape(shape)
# du / dt  [e / s] = [kg m2 / s3]
rate = rate.reshape(shape) * units.UnitMass_in_cgs * units.UnitLength_in_cgs**2 / units.UnitTime_in_cgs**3
cs = cs.reshape(shape) * units.UnitLength_in_cgs / units.UnitTime_in_cgs

L_cool = -rate * cs / astropy.units.kpc.to("cm")

plt.title("Cooling Length")
plt.contourf(rho, T, np.log(L_cool))
plt.xlabel("Density [$cm^{-3}$]")
plt.ylabel("Temperature [K]")
plt.colorbar()

ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()
