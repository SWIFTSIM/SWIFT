# Script for plotting energy evolution of uniform box of gas with single star in the 
# centre when running with stochastic energy injection. It also checks that the change
# in total energy of the gas particles is within a specified tolerance from what is 
# expected based on the mass of the star particle (Note that this tolerance could be
# somewhat high because of Poisson noise and the relatively small number of injection
# events)

import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py
import os.path
import numpy as np
import glob

# Number of snapshots and elements
newest_snap_name = max(glob.glob('stellar_evolution_*.hdf5'), key=os.path.getctime)
n_snapshots = int(newest_snap_name.replace('stellar_evolution_','').replace('.hdf5','')) + 1
n_elements = 9

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 9,
'legend.fontsize': 9,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.3,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.18,
'figure.subplot.top'     : 0.92,
'figure.subplot.wspace'  : 0.21,
'figure.subplot.hspace'  : 0.19,
'lines.markersize' : 6,
'lines.linewidth' : 2.,
'text.latex.unicode': True
}

rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})

# Read the simulation data
sim = h5py.File("stellar_evolution_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"][0]
kernel = sim["/HydroScheme"].attrs["Kernel function"][0]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = sim["/HydroScheme"].attrs["Kernel eta"][0]
alpha = sim["/HydroScheme"].attrs["Alpha viscosity"][0]
H_mass_fraction = sim["/HydroScheme"].attrs["Hydrogen mass fraction"][0]
H_transition_temp = sim["/HydroScheme"].attrs["Hydrogen ionization transition temperature"][0]
T_initial = sim["/HydroScheme"].attrs["Initial temperature"][0]
T_minimal = sim["/HydroScheme"].attrs["Minimal temperature"][0]
git = sim["Code"].attrs["Git Revision"]
star_initial_mass = sim["/PartType4/Masses"][0]

# Cosmological parameters
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
gas_gamma = sim["/HydroScheme"].attrs["Adiabatic index"][0]

# Units
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
unit_vel_in_cgs = unit_length_in_cgs / unit_time_in_cgs
unit_energy_in_cgs = unit_mass_in_cgs * unit_vel_in_cgs * unit_vel_in_cgs
unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs

# Calculate solar mass in internal units
const_solar_mass = 1.98848e33 / unit_mass_in_cgs

# Define Gyr
Gyr_in_cgs = 1e9 * 365 * 24 * 3600.

# Find out how many particles (gas and star) we have
n_parts = sim["/Header"].attrs["NumPart_Total"][0]
n_sparts = sim["/Header"].attrs["NumPart_Total"][4]

# Declare arrays for data
masses = zeros((n_parts,n_snapshots))
star_masses = zeros((n_sparts,n_snapshots))
internal_energy = zeros((n_parts,n_snapshots))
velocity_parts = zeros((n_parts,3,n_snapshots))
time = zeros(n_snapshots)

# Read fields we are checking from snapshots
#for i in [0,n_snapshots-1]:
for i in range(n_snapshots):
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	print('reading snapshot '+str(i))
	masses[:,i] = sim["/PartType0/Masses"]
	internal_energy[:,i] = sim["/PartType0/InternalEnergy"]
	velocity_parts[:,:,i] = sim["/PartType0/Velocities"]
	time[i] = sim["/Header"].attrs["Time"][0]

# Check that the total amount of enrichment is as expected.
# Define tolerance. Note, relatively high value used due to
# Poisson noise associated with stochastic energy injection.
eps = 0.15

# Stochastic heating
vel2 = zeros((n_parts,n_snapshots))
vel2[:,:] = velocity_parts[:,0,:]*velocity_parts[:,0,:] + velocity_parts[:,1,:]*velocity_parts[:,1,:] + velocity_parts[:,2,:]*velocity_parts[:,2,:]
total_kinetic_energy_cgs = np.sum(np.multiply(vel2,masses)*0.5,axis = 0) * unit_energy_in_cgs
total_energy_cgs = np.sum(np.multiply(internal_energy,masses),axis = 0) * unit_energy_in_cgs
total_energy_released_cgs = total_energy_cgs[n_snapshots-1] - total_energy_cgs[0] + total_kinetic_energy_cgs[n_snapshots-1] - total_kinetic_energy_cgs[0]

# Calculate energy released
energy_per_sn = 1.0e51 / unit_energy_in_cgs
SNIa_efficiency = 2.e-3
SNIa_timescale_Gyr = 2.0
expected_energy_released_cgs = np.zeros(n_snapshots)
for i in range(n_snapshots):
	age_Gyr = time[i] * unit_time_in_cgs / Gyr_in_cgs
	total_sn = SNIa_efficiency * (SNIa_timescale_Gyr * Gyr_in_cgs / unit_time_in_cgs * (1 - exp(-age_Gyr/SNIa_timescale_Gyr))) * star_initial_mass * const_solar_mass
	expected_energy_released_cgs[i] = total_sn * energy_per_sn * unit_energy_in_cgs

# Did we get it right?
if abs(total_energy_released_cgs - expected_energy_released_cgs[n_snapshots-1])/expected_energy_released_cgs[n_snapshots-1] < eps:
	print("total stochastic energy release consistent with expectation. total stochastic energy release "+str(total_energy_released_cgs)+" expected "+ str(expected_energy_released_cgs[n_snapshots-1]) + " initial total internal energy "+ str(total_energy_cgs[0] + total_kinetic_energy_cgs[0]))
else:
	print("total stochastic energy release "+str(total_energy_released_cgs)+" expected "+ str(expected_energy_released_cgs[n_snapshots-1]) + " initial total internal energy "+ str(total_energy_cgs[0] + total_kinetic_energy_cgs[0]) + " energy change fraction of total " + str(total_energy_released_cgs/(total_energy_cgs[0]+total_kinetic_energy_cgs[0])))

# Plot the energy evolution
figure()
subplot(111)
plot(time*unit_time_in_cgs/Gyr_in_cgs, total_energy_cgs + total_kinetic_energy_cgs - total_energy_cgs[0] - total_kinetic_energy_cgs[0],color='k', linewidth=0.5, label="SWIFT")
plot(time*unit_time_in_cgs/Gyr_in_cgs, expected_energy_released_cgs,color = 'r', linewidth=0.5, label="expected")
xlabel("Time (Gyr)")
ylabel("Total energy (erg)")
legend()
savefig("continuous_energy_evolution.png", dpi=200)

