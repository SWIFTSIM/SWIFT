import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import sys
import glob

# Get the total number of snapshots
file_list = glob.glob("CoolingHalo_*")
n_snaps = len(file_list)

# some constants
OMEGA = 0.3  # Cosmological matter fraction at z = 0
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
h = 0.67777  # hubble parameter
gamma = 5.0 / 3.0
eta = 1.2349
H_0_cgs = 100.0 * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

# read some header/parameter information from the first snapshot

filename = "CoolingHalo_0000.hdf5"
f = h5.File(filename, "r")
params = f["Parameters"]
unit_mass_cgs = float(params.attrs["InternalUnitSystem:UnitMass_in_cgs"])
unit_length_cgs = float(params.attrs["InternalUnitSystem:UnitLength_in_cgs"])
unit_velocity_cgs = float(params.attrs["InternalUnitSystem:UnitVelocity_in_cgs"])
unit_time_cgs = unit_length_cgs / unit_velocity_cgs
v_c = float(params.attrs["IsothermalPotential:vrot"])
v_c_cgs = v_c * unit_velocity_cgs
header = f["Header"]
N = header.attrs["NumPart_Total"][0]
box_centre = np.array(header.attrs["BoxSize"])

# calculate r_vir and M_vir from v_c
r_vir_cgs = v_c_cgs / (10.0 * H_0_cgs * np.sqrt(OMEGA))
M_vir_cgs = r_vir_cgs * v_c_cgs ** 2 / CONST_G_CGS

potential_energy_array = np.zeros(n_snaps)
internal_energy_array = np.zeros(n_snaps)
kinetic_energy_array = np.zeros(n_snaps)
time_array_cgs = np.zeros(n_snaps)

s_potential_energy_array = np.zeros(n_snaps)
s_internal_energy_array = np.zeros(n_snaps)
s_kinetic_energy_array = np.zeros(n_snaps)
s_rad_energy_array = np.zeros(n_snaps)

# get the radiated energy
statistics = np.loadtxt("statistics.txt")
stats_kinetic_energy_array = statistics[:, 13]
stats_potential_energy_array = statistics[:, 14]
stats_internal_energy_array = statistics[:, 15]
stats_rad_energy_array = statistics[:, 16]
stats_times_array = statistics[:, 1]

for i in range(n_snaps):

    filename = "CoolingHalo_%04d.hdf5" % i
    f = h5.File(filename, "r")
    coords_dset = f["PartType0/Coordinates"]
    coords = np.array(coords_dset)
    # translate coords by centre of box
    header = f["Header"]
    snap_time = header.attrs["Time"]
    snap_time_cgs = snap_time * unit_time_cgs
    time_array_cgs[i] = snap_time_cgs

    coords[:, 0] -= box_centre[0] / 2.0
    coords[:, 1] -= box_centre[1] / 2.0
    coords[:, 2] -= box_centre[2] / 2.0
    radius = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2 + coords[:, 2] ** 2)
    radius_cgs = radius * unit_length_cgs
    radius_over_virial_radius = radius_cgs / r_vir_cgs

    r = radius_over_virial_radius
    total_potential_energy = np.sum(v_c ** 2 * np.log(r))
    potential_energy_array[i] = total_potential_energy

    vels_dset = f["PartType0/Velocities"]
    vels = np.array(vels_dset)
    speed_squared = vels[:, 0] ** 2 + vels[:, 1] ** 2 + vels[:, 2] ** 2
    total_kinetic_energy = 0.5 * np.sum(speed_squared)
    kinetic_energy_array[i] = total_kinetic_energy

    u_dset = f["PartType0/InternalEnergies"]
    u = np.array(u_dset)
    total_internal_energy = np.sum(u)
    internal_energy_array[i] = total_internal_energy

    # find radiated energy index based on current time in statistics.txt
    time_ind = np.argmin(np.abs(stats_times_array - snap_time))
    s_kinetic_energy_array[i] = stats_kinetic_energy_array[time_ind]
    s_potential_energy_array[i] = stats_potential_energy_array[time_ind]
    s_internal_energy_array[i] = stats_internal_energy_array[time_ind]
    s_rad_energy_array[i] = stats_rad_energy_array[time_ind]


# put energies in units of v_c^2 and rescale by number of particles

pe = potential_energy_array / (N * v_c ** 2)
ke = kinetic_energy_array / (N * v_c ** 2)
ie = internal_energy_array / (N * v_c ** 2)
te = pe + ke + ie  # + re

print(pe)
print(ke)
print(ie)
#  print(re)
print(te)

spe = s_potential_energy_array / (N * v_c ** 2)
ske = s_kinetic_energy_array / (N * v_c ** 2)
sie = s_internal_energy_array / (N * v_c ** 2)
sre = s_rad_energy_array / (N * v_c ** 2)
ste = spe + ske + sie + sre

dyn_time_cgs = r_vir_cgs / v_c_cgs
time_array = time_array_cgs / dyn_time_cgs

fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax1.plot(time_array, ke, label="Kinetic Energy")
ax1.plot(time_array, pe, label="Potential Energy")
ax1.plot(time_array, ie, label="Internal Energy")
ax1.plot(time_array, te, label="Total Energy")

ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(time_array, ske, label="Kinetic Energy")
ax2.plot(time_array, spe, label="Potential Energy")
ax2.plot(time_array, sie, label="Internal Energy")
ax2.plot(time_array, sre, label="Radiated Energy")
ax2.plot(time_array, ste, label="Total Energy")

for ax in [ax1, ax2]:
    ax.legend(loc="lower left")
    ax.set_xlabel(r"$t / t_{dyn}$")
    ax.set_ylabel(r"$E / v_c^2$")
ax1.set_title("From particles")
ax2.set_title("From statistics")

fig.suptitle(
    r"$%d \, \, \mathrm{particles} \,,\, v_c = %.1f \, \mathrm{km / s}$" % (N, v_c)
)
# plt.ylim((-4,2))
# plot_filename = "density_profile_%03d.png" %i
plt.show()
