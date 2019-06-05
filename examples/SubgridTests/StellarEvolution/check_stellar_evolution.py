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
energy_per_SNe = 1.0e51/unit_energy_in_cgs

# Check that the total amount of enrichment is as expected.
# Define tolerance
eps = 0.01

# Total mass
total_part_mass = np.sum(masses,axis = 0)
if abs((total_part_mass[n_snapshots-1] - total_part_mass[0])/total_part_mass[0] - ejected_mass/total_part_mass[0])*total_part_mass[0]/ejected_mass < eps:
	print("total mass released consistent with expectation")
else:
	print("mass increase "+str(total_part_mass[n_snapshots-1]/total_part_mass[0])+" expected "+ str(1.0+ejected_mass/total_part_mass[0]))

# Check that mass is conserved (i.e. total star mass decreases by same amount as total gas mass increases)
total_spart_mass = np.sum(star_masses,axis = 0)
if abs((total_part_mass[n_snapshots-1] + total_spart_mass[n_snapshots-1]) / (total_part_mass[0] + total_spart_mass[0]) - 1.0) < eps**3:
	print("total mass conserved")
else:
	print("initial part, spart mass " + str(total_part_mass[0]) + " " + str(total_spart_mass[0]) + " final mass " + str(total_part_mass[n_snapshots-1]) + " " + str(total_spart_mass[n_snapshots-1]))

# Total metal mass from AGB
total_metal_mass_AGB = np.sum(np.multiply(metal_mass_frac_from_AGB,masses),axis = 0)
expected_metal_mass_AGB = ejecta_factor*ejected_mass
if abs(total_metal_mass_AGB[n_snapshots-1] - expected_metal_mass_AGB)/expected_metal_mass_AGB < eps:
	print("total AGB metal mass released consistent with expectation")
else:
	print("total AGB metal mass "+str(total_metal_mass_AGB[n_snapshots-1])+" expected "+ str(expected_metal_mass_AGB))

# Total mass from AGB
total_AGB_mass = np.sum(mass_from_AGB,axis = 0)
expected_AGB_mass = ejecta_factor*ejected_mass
if abs(total_AGB_mass[n_snapshots-1] - expected_AGB_mass)/expected_AGB_mass < eps:
	print("total AGB mass released consistent with expectation")
else:
	print("total AGB mass "+str(total_AGB_mass[n_snapshots-1])+" expected "+ str(expected_AGB_mass))

# Total metal mass from SNII
total_metal_mass_SNII = np.sum(np.multiply(metal_mass_frac_from_SNII,masses),axis = 0)
expected_metal_mass_SNII = ejecta_factor*ejected_mass
if abs(total_metal_mass_SNII[n_snapshots-1] - expected_metal_mass_SNII)/expected_metal_mass_SNII < eps:
	print("total SNII metal mass released consistent with expectation")
else:
	print("total SNII metal mass "+str(total_metal_mass_SNII[n_snapshots-1])+" expected "+ str(expected_metal_mass_SNII))

# Total mass from SNII
total_SNII_mass = np.sum(mass_from_SNII,axis = 0)
expected_SNII_mass = ejecta_factor*ejected_mass
if abs(total_SNII_mass[n_snapshots-1] - expected_SNII_mass)/expected_SNII_mass < eps:
	print("total SNII mass released consistent with expectation")
else:
	print("total SNII mass "+str(total_SNII_mass[n_snapshots-1])+" expected "+ str(expected_SNII_mass))

# Total metal mass from SNIa
total_metal_mass_SNIa = np.sum(np.multiply(metal_mass_frac_from_SNIa,masses),axis = 0)
expected_metal_mass_SNIa = ejecta_factor*ejected_mass
if abs(total_metal_mass_SNIa[n_snapshots-1] - expected_metal_mass_SNIa)/expected_metal_mass_SNIa < eps:
	print("total SNIa metal mass released consistent with expectation")
else:
	print("total SNIa metal mass "+str(total_metal_mass_SNIa[n_snapshots-1])+" expected "+ str(expected_metal_mass_SNIa))

# Total iron mass from SNIa
total_iron_mass_SNIa = np.sum(np.multiply(iron_mass_frac_from_SNIa,masses),axis = 0)
expected_iron_mass_SNIa = ejecta_factor*ejected_mass
if abs(total_iron_mass_SNIa[n_snapshots-1] - expected_iron_mass_SNIa)/expected_iron_mass_SNIa < eps:
	print("total SNIa iron mass released consistent with expectation")
else:
	print("total SNIa iron mass "+str(total_iron_mass_SNIa[n_snapshots-1])+" expected "+ str(expected_iron_mass_SNIa))

# Total mass from SNIa
total_SNIa_mass = np.sum(mass_from_SNIa,axis = 0)
expected_SNIa_mass = ejecta_factor*ejected_mass
if abs(total_SNIa_mass[n_snapshots-1] - expected_SNIa_mass)/expected_SNIa_mass < eps:
	print("total SNIa mass released consistent with expectation")
else:
	print("total SNIa mass "+str(total_SNIa_mass[n_snapshots-1])+" expected "+ str(expected_SNIa_mass))

# Total metal mass
total_metal_mass = np.sum(np.multiply(metallicity,masses),axis = 0)
expected_metal_mass = ejecta_factor_metallicity*ejected_mass
if abs(total_metal_mass[n_snapshots-1] - expected_metal_mass)/expected_metal_mass < eps:
	print("total metal mass released consistent with expectation")
else:
	print("total metal mass "+str(total_metal_mass[n_snapshots-1])+" expected "+ str(expected_metal_mass))

# Total mass for each element
expected_element_mass = ejecta_factor_abundances*ejected_mass
for i in range(n_elements):
	total_element_mass = np.sum(np.multiply(abundances[:,i,:],masses),axis = 0)
	if abs(total_element_mass[n_snapshots-1] - expected_element_mass)/expected_element_mass < eps:
		print("total element mass released consistent with expectation for element "+str(i))
	else:
		print("total element mass "+str(total_element_mass[n_snapshots-1])+" expected "+ str(expected_element_mass) + " for element "+ str(i))
