import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py
import os.path
import numpy as np

def sort_by_id(id,x):
	indices = np.argsort(id)
	x_sorted = zeros(len(x))
	for i in range(len(x)):
		x_sorted[i] = x[indices[i]]
	return x_sorted

def distance(p,s):
	dist2 = (p[0] - s[0]) * (p[0] - s[0]) + (p[1] - s[1]) * (p[1] - s[1]) +(p[2] - s[2]) * (p[2] - s[2])
	return sqrt(dist2)

n_snapshots = 10

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 9,
'legend.fontsize': 9,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.14,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.12,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
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

# Cosmological parameters
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
gas_gamma = sim["/HydroScheme"].attrs["Adiabatic index"][0]

unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]

unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs

n_parts = sim["/Header"].attrs["NumPart_Total"][0]
n_sparts = sim["/Header"].attrs["NumPart_Total"][4]

metallicity = zeros((n_parts,n_snapshots))
coord_parts = zeros((n_parts,3))
coord_sparts = zeros(3)
id = zeros((n_parts,n_snapshots))
time = zeros(n_snapshots)

for i in range(n_parts):
	id[:,0] = sim["/PartType0/ParticleIDs"]
	coord_parts_unsorted = sim["/PartType0/Coordinates"]
	coord_parts = sort_by_id(id[:,0], coord_parts_unsorted)

if n_sparts == 1:
	coord_sparts = sim["/PartType4/Coordinates"]
	smoothing_length = sim["/PartType4/SmoothingLength"]
else:
	print("too many sparts for now")

for i in range(n_snapshots):
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	id[:,i] = sim["/PartType0/ParticleIDs"]
	metallicity_unsorted = sim["/PartType0/Metallicity"]
	metallicity[:,i] = sort_by_id(id[:,i], metallicity_unsorted)
	time[i] = sim["/Header"].attrs["Time"][0]

figure()
subplot(111)

for i in range(100):
	if (distance(coord_parts[i,:],coord_sparts) < smoothing_length):
		plot(time,metallicity[i,:],color='k',alpha=0.1)

xlabel("time")
ylabel("metallicity")
savefig("stellar_evolution.png", dpi=200)
