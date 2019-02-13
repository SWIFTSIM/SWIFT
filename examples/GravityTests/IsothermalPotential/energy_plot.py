import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.145,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.11,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


import numpy as np
import h5py as h5
import sys

# File containing the total energy
stats_filename = "./energy.txt"

# First snapshot
snap_filename = "Isothermal_0000.hdf5"
f = h5.File(snap_filename,'r')

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]

# Read the header
header = f["Header"]
box_size = float(header.attrs["BoxSize"][0])

# Read the properties of the potential
parameters = f["Parameters"]
R200 = 100 
Vrot = float(parameters.attrs["IsothermalPotential:vrot"])
centre = [box_size/2, box_size/2, box_size/2]
f.close()

# Read the statistics summary
file_energy = np.loadtxt("energy.txt")
time_stats = file_energy[:,0]
E_kin_stats = file_energy[:,3]
E_pot_stats = file_energy[:,5]
E_tot_stats = E_kin_stats + E_pot_stats

# Read the snapshots
time_snap = np.zeros(402)
E_kin_snap = np.zeros(402)
E_pot_snap = np.zeros(402)
E_tot_snap = np.zeros(402)
Lz_snap = np.zeros(402)

# Read all the particles from the snapshots
for i in range(402):
    snap_filename = "Isothermal_%0.4d.hdf5"%i
    f = h5.File(snap_filename,'r')

    pos_x = f["PartType1/Coordinates"][:,0]
    pos_y = f["PartType1/Coordinates"][:,1]
    pos_z = f["PartType1/Coordinates"][:,2]
    vel_x = f["PartType1/Velocities"][:,0]
    vel_y = f["PartType1/Velocities"][:,1]
    vel_z = f["PartType1/Velocities"][:,2]
    mass = f["/PartType1/Masses"][:]
    
    r = np.sqrt((pos_x[:] - centre[0])**2 + (pos_y[:] - centre[1])**2 + (pos_z[:] - centre[2])**2)
    Lz = (pos_x[:] - centre[0]) * vel_y[:] - (pos_y[:] - centre[1]) * vel_x[:]

    time_snap[i] = f["Header"].attrs["Time"]
    E_kin_snap[i] = np.sum(0.5 * mass * (vel_x[:]**2 + vel_y[:]**2 + vel_z[:]**2))
    E_pot_snap[i] = np.sum(mass * Vrot**2 *  log(r))
    E_tot_snap[i] = E_kin_snap[i] + E_pot_snap[i]
    Lz_snap[i] = np.sum(Lz)

# Plot energy evolution
figure()
plot(time_stats, E_kin_stats, "r-", lw=0.5, label="Kinetic energy")
plot(time_stats, E_pot_stats, "g-", lw=0.5, label="Potential energy")
plot(time_stats, E_tot_stats, "k-", lw=0.5, label="Total energy")

plot(time_snap[::10], E_kin_snap[::10], "rD", lw=0.5, ms=2)
plot(time_snap[::10], E_pot_snap[::10], "gD", lw=0.5, ms=2)
plot(time_snap[::10], E_tot_snap[::10], "kD", lw=0.5, ms=2)

legend(loc="center right", fontsize=8, frameon=False, handlelength=3, ncol=1)
xlabel("${\\rm{Time}}$", labelpad=0)
ylabel("${\\rm{Energy}}$",labelpad=0)
xlim(0, 8)

savefig("energy.png", dpi=200)

# Plot angular momentum evolution
figure()
plot(time_snap, Lz_snap, "k-", lw=0.5, ms=2)

xlabel("${\\rm{Time}}$", labelpad=0)
ylabel("${\\rm{Angular~momentum}}$",labelpad=0)
xlim(0, 8)

savefig("angular_momentum.png", dpi=200)


