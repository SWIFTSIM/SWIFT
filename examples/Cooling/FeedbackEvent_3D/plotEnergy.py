"""
Plots the energy from the statistics.txt file for this simulation.
"""

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import load

from unyt import Gyr, erg, mh, kb

from makeIC import (
    gamma,
    initial_temperature,
    inject_temperature,
    mu,
    particle_mass,
)

try:
    plt.style.use("mnras_durham")
except:
    pass


# Snapshot for grabbing the units.
snapshot = load("feedback_0000.hdf5")
units = snapshot.metadata.units
energy_units = units.mass * units.length ** 2 / (units.time ** 2)

data = np.loadtxt("statistics.txt").T

#  Assign correct units to each

time = data[1] * units.time
mass = data[5] * units.mass
kinetic_energy = data[13] * energy_units
thermal_energy = data[14] * energy_units
potential_energy = data[15] * energy_units
radiative_cool = data[16] * energy_units
total_energy = kinetic_energy + thermal_energy + potential_energy

# Now we have to figure out how much energy we actually 'injected'
background_internal_energy = (
    (1.0 / (mu * mh)) * (kb / (gamma - 1.0)) * initial_temperature
)

heated_internal_energy = (1.0 / (mu * mh)) * (kb / (gamma - 1.0)) * inject_temperature

injected_energy = (heated_internal_energy - background_internal_energy) * particle_mass

#  Also want to remove the 'background' energy
n_parts = snapshot.metadata.n_gas
total_background_energy = background_internal_energy * n_parts * particle_mass

# Now we can plot

fig, ax = plt.subplots()

ax.plot(time.to(Gyr), (kinetic_energy).to(erg), label="Kinetic")

ax.plot(
    time.to(Gyr), (thermal_energy - total_background_energy).to(erg), label="Thermal"
)

ax.plot(time.to(Gyr), (radiative_cool).to(erg), label="Lost to cooling")

ax.set_xlim(0, 0.05 * Gyr)

ax.set_xlabel("Time [Gyr]")
ax.set_ylabel("Energy [erg]")

ax.legend()

fig.tight_layout()
fig.savefig("Energy.pdf")
