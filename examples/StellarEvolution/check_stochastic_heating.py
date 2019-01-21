import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py
import os.path
import numpy as np

# Number of snapshots and elements
n_snapshots = 11
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

# Find out how many particles (gas and star) we have
n_parts = sim["/Header"].attrs["NumPart_Total"][0]
n_sparts = sim["/Header"].attrs["NumPart_Total"][4]

# Declare arrays for data
masses = zeros((n_parts,n_snapshots))
star_masses = zeros((n_sparts,n_snapshots))
mass_from_AGB = zeros((n_parts,n_snapshots))
metal_mass_frac_from_AGB = zeros((n_parts,n_snapshots))
mass_from_SNII = zeros((n_parts,n_snapshots))
metal_mass_frac_from_SNII = zeros((n_parts,n_snapshots))
mass_from_SNIa = zeros((n_parts,n_snapshots))
metal_mass_frac_from_SNIa = zeros((n_parts,n_snapshots))
iron_mass_frac_from_SNIa = zeros((n_parts,n_snapshots))
metallicity = zeros((n_parts,n_snapshots))
abundances = zeros((n_parts,n_elements,n_snapshots))
internal_energy = zeros((n_parts,n_snapshots))
coord_parts = zeros((n_parts,3))
velocity_parts = zeros((n_parts,3,n_snapshots))
coord_sparts = zeros(3)
time = zeros(n_snapshots)

# Read fields we are checking from snapshots
for i in range(n_snapshots):
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	print('reading snapshot '+str(i))
	abundances[:,:,i] = sim["/PartType0/ElementAbundance"]
	metallicity[:,i] = sim["/PartType0/Metallicity"]
	masses[:,i] = sim["/PartType0/Masses"]
	star_masses[:,i] = sim["/PartType4/Masses"]
	mass_from_AGB[:,i] = sim["/PartType0/TotalMassFromAGB"]
	metal_mass_frac_from_AGB[:,i] = sim["/PartType0/MetalMassFracFromAGB"]
	mass_from_SNII[:,i] = sim["/PartType0/TotalMassFromSNII"]
	metal_mass_frac_from_SNII[:,i] = sim["/PartType0/MetalMassFracFromSNII"]
	mass_from_SNIa[:,i] = sim["/PartType0/TotalMassFromSNIa"]
	metal_mass_frac_from_SNIa[:,i] = sim["/PartType0/MetalMassFracFromSNIa"]
	iron_mass_frac_from_SNIa[:,i] = sim["/PartType0/IronMassFracFromSNIa"]
	internal_energy[:,i] = sim["/PartType0/InternalEnergy"]
	velocity_parts[:,:,i] = sim["/PartType0/Velocities"]
	time[i] = sim["/Header"].attrs["Time"][0]

# Define ejecta factor
ejecta_factor = 1.0e-2
ejecta_factor_metallicity = 1.0 - 2.0/n_elements
ejecta_factor_abundances = 1.0/n_elements
ejected_mass = star_initial_mass
SNIa_rate = 1.0e9
energy_per_SNe = 1.0e51/unit_energy_in_cgs

# Check that the total amount of enrichment is as expected.
# Define tolerance
eps = 0.01

# print mean, max internal energy of particles
for i in range(n_snapshots):
	mean_u_cgs = np.mean(internal_energy[:,i]*unit_vel_in_cgs*unit_vel_in_cgs)
	max_u_cgs = np.max(internal_energy[:,i]*unit_vel_in_cgs*unit_vel_in_cgs)
	print("snapshot "+ str(i) + " mean max u cgs " + str(mean_u_cgs) + " " + str(max_u_cgs))


# Continuous heating
#vel2 = zeros((n_parts,n_snapshots))
#vel2[:,:] = velocity_parts[:,0,:]*velocity_parts[:,0,:] + velocity_parts[:,1,:]*velocity_parts[:,1,:] + velocity_parts[:,2,:]*velocity_parts[:,2,:]
#total_kinetic_energy = np.sum(np.multiply(vel2,masses)*0.5,axis = 0)
#total_energy = np.sum(np.multiply(internal_energy,masses),axis = 0)
#total_energy_released = total_energy[n_snapshots-1] - total_energy[0] + total_kinetic_energy[n_snapshots-1] - total_kinetic_energy[0]
#expected_energy_released = SNIa_rate*star_initial_mass*time[n_snapshots-1]*energy_per_SNe
#if abs(total_energy_released - expected_energy_released)/expected_energy_released < eps:
#	print("total continuous energy release consistent with expectation")
#else:
#	print("total continuous energy release "+str(total_energy_released)+" expected "+ str(expected_energy_released) + " initial total internal energy "+ str(total_energy[0] + total_kinetic_energy[0]) + " energy change fraction of total " + str(total_energy_released/(total_energy[0]+total_kinetic_energy[0])))

# Stochastic heating
#vel2 = zeros((n_parts,n_snapshots))
#vel2[:,:] = velocity_parts[:,0,:]*velocity_parts[:,0,:] + velocity_parts[:,1,:]*velocity_parts[:,1,:] + velocity_parts[:,2,:]*velocity_parts[:,2,:]
#total_kinetic_energy = np.sum(np.multiply(vel2,masses)*0.5,axis = 0)
#total_energy = np.sum(np.multiply(internal_energy,masses),axis = 0)
#total_energy_released = total_energy[n_snapshots-1] - total_energy[0] + total_kinetic_energy[n_snapshots-1] - total_kinetic_energy[0]
#
## put in variable names to make this clearer...
#heating_probability = 2.39802e-01 * 9.76562e+01 / (3.16228e+07 * 2e-4)
#du = 3.16228e+07 * 2.0e-2
#print(heating_probability,du)
#
#expected_energy_released = heating_probability * du * np.sum(star_masses,axis=0)[0]
#if abs(total_energy_released - expected_energy_released)/expected_energy_released < eps:
#	print("total stochastic energy release consistent with expectation")
#else:
#	print("total stochastic energy release "+str(total_energy_released)+" expected "+ str(expected_energy_released) + " initial total internal energy "+ str(total_energy[0] + total_kinetic_energy[0]) + " energy change fraction of total " + str(total_energy_released/(total_energy[0]+total_kinetic_energy[0])))
