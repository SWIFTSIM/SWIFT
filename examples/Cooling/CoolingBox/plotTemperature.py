from h5py import File
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Plot parameters
params = {
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'font.size': 12,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': (3.15, 3.15),
    'figure.subplot.left': 0.145,
    'figure.subplot.right': 0.99,
    'figure.subplot.bottom': 0.11,
    'figure.subplot.top': 0.99,
    'figure.subplot.wspace': 0.15,
    'figure.subplot.hspace': 0.12,
    'lines.markersize': 6,
    'lines.linewidth': 3.,
}
plt.rcParams.update(params)


# Some constants in cgs units
k_b_cgs = 1.38e-16  # boltzmann
m_h_cgs = 1.67e-24  # proton mass
# need to be changed in makeIC.py too
h_frac = 0.76
mu = 4. / (1. + 3. * h_frac)


# File containing the total energy
stats_filename = "./energy.txt"

# First snapshot
snap_filename = "coolingBox_0000.hdf5"

# Read the initial state of the gas
f = File(snap_filename, 'r')

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]

# Read the adiabatic index
gamma = float(f["HydroScheme"].attrs["Adiabatic index"])


def Temperature(u):
    """ Compute the temperature from the internal energy. """
    u *= (unit_length / unit_time)**2
    return u * (gamma - 1.) * m_h_cgs / (mu * k_b_cgs)


# Read energy and time arrays
array = np.genfromtxt(stats_filename, skip_header=1)
time = array[:, 0] * unit_time
total_mass = array[:, 1]
total_energy = array[:, 2]
kinetic_energy = array[:, 3]
internal_energy = array[:, 4]
radiated_energy = array[:, 8]
initial_energy = total_energy[0]

# Conversions to cgs
total_energy_cgs = total_energy / total_mass[0]
total_energy_cgs = Temperature(total_energy_cgs)

kinetic_energy_cgs = kinetic_energy / total_mass[0]
kinetic_energy_cgs = Temperature(kinetic_energy_cgs)

internal_energy_cgs = internal_energy / total_mass[0]
internal_energy_cgs = Temperature(internal_energy_cgs)

radiated_energy_cgs = radiated_energy / total_mass[0]
radiated_energy_cgs = Temperature(radiated_energy_cgs)

# Read snapshots
temp_snap = np.zeros(25)
time_snap_cgs = np.zeros(25)
for i in range(25):
    snap = File("coolingBox_%0.4d.hdf5" % i, 'r')
    u = snap["/PartType0/InternalEnergy"][:] * snap["/PartType0/Masses"][:]
    u = sum(u) / total_mass[0]
    temp_snap[i] = Temperature(u)
    time_snap_cgs[i] = snap["/Header"].attrs["Time"] * unit_time


plt.figure()

Myr_in_yr = 3.15e13
plt.plot(time, total_energy_cgs, 'r-', lw=1.6, label="Gas total temperature")
plt.plot(time_snap_cgs, temp_snap, 'rD', ms=3)
plt.plot(time, radiated_energy_cgs, 'g-', lw=1.6, label="Radiated temperature")
plt.plot(time, total_energy_cgs + radiated_energy_cgs, 'b-',
         lw=0.6, label="Gas total + radiated")

plt.legend(loc="upper right", fontsize=8, frameon=False,
           handlelength=3, ncol=1)
plt.xlabel("${\\rm{Time~[Myr]}}$", labelpad=0)
plt.ylabel("${\\rm{Temperature~[K]}}$")

plt.savefig("temperature.png", dpi=200)
