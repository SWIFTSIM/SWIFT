"""
Plots the energy from the energy.txt file for this simulation.
"""

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import load

from unyt import Gyr, erg, mh, kb, Myr
from scipy.interpolate import interp1d

from makeIC import gamma, initial_density, initial_temperature, inject_temperature, mu, particle_mass

try:
    plt.style.use("mnras_durham")
except:
    pass

time_to_plot = 25 * Myr
diffusion_parameters = [0.1 * x for x in range(11)]
plot_directory_name = "default_diffmax"

kinetic_energy_at_time = []
thermal_energy_at_time = []
radiative_energy_at_time = []

for diffusion in diffusion_parameters:
    directory_name = f"{plot_directory_name}_{diffusion:1.1f}"

# Snapshot for grabbing the units.
    snapshot = load(f"{directory_name}/feedback_0000.hdf5")
    units = snapshot.metadata.units
    energy_units = units.mass * units.length ** 2 / (units.time ** 2)

    data = np.loadtxt(f"{directory_name}/energy.txt").T

# Assign correct units to each

    time = data[0] * units.time
    mass = data[1] * units.mass
    total_energy = data[2] * energy_units
    kinetic_energy = data[3] * energy_units
    thermal_energy = data[4] * energy_units
    radiative_cool = data[8] * energy_units

# Now we have to figure out how much energy we actually 'injected'
    background_internal_energy = (
        (1.0 / (mu * mh)) * (kb / (gamma - 1.0)) * initial_temperature
    )

    heated_internal_energy = (1.0 / (mu * mh)) * (kb / (gamma - 1.0)) * inject_temperature

    injected_energy = (heated_internal_energy - background_internal_energy) * particle_mass

# Also want to remove the 'background' energy
    n_parts = snapshot.metadata.n_gas
    total_background_energy = background_internal_energy * n_parts * particle_mass
    
    kinetic_energy_interpolated = interp1d(
        time.to(Myr),
        kinetic_energy.to(erg)
    )

    thermal_energy_interpolated = interp1d(
        time.to(Myr),
        (thermal_energy - total_background_energy).to(erg)
    )

    radiative_cool_interpolated = interp1d(
        time.to(Myr),
        radiative_cool.to(erg)
    )

    kinetic_energy_at_time.append(kinetic_energy_interpolated(time_to_plot.to(Myr)))
    thermal_energy_at_time.append(thermal_energy_interpolated(time_to_plot.to(Myr)))
    radiative_energy_at_time.append(radiative_cool_interpolated(time_to_plot.to(Myr)))


# Now we can plot


fig, ax = plt.subplots()

ax.plot(
    diffusion_parameters,
    kinetic_energy_at_time,
    label="Kinetic"
)

ax.plot(
    diffusion_parameters,
    thermal_energy_at_time,
    label="Thermal"
)

ax.plot(
    diffusion_parameters,
    radiative_energy_at_time,
    label="Lost to cooling"
)

ax.set_xlim(0, 1.0)

ax.set_xlabel(r"Diffusion $\alpha_{\rm max}$")
ax.set_ylabel(f"Energy in component at $t=${time_to_plot} [erg]")

ax.legend()

fig.tight_layout()
fig.savefig("EnergyFuncDiff.pdf")
