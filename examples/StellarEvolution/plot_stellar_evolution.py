import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py
import os.path
import numpy as np

# Function to sort the coordinates of particles based on their ids
# id: array of particle ids
# x: array of coordinates
def sort_coords_by_id(id,x):
	indices = np.argsort(id)
	x_sorted = zeros(x.shape)
	for i in range(len(x)):
		x_sorted[i,:] = x[indices[i],:]
	return x_sorted

# Function to sort particle data based on their ids
# id: array of particle ids
# x: array of particle data
def sort_by_id(id,x):
	indices = np.argsort(id)
	x_sorted = zeros(x.shape)
	for i in range(len(x)):
		x_sorted[i] = x[indices[i]]
	return x_sorted

# Function to determine distance between part and spart
# p: part coordinates
# s: spart coordinates
def distance(p,s):
	dist2 = (p[0] - s[0]) * (p[0] - s[0]) + (p[1] - s[1]) * (p[1] - s[1]) +(p[2] - s[2]) * (p[2] - s[2])
	return sqrt(dist2)

# Number of snapshots we have
n_snapshots = 11

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

# Find out how many particles (gas and star) we have
n_parts = sim["/Header"].attrs["NumPart_Total"][0]
n_sparts = sim["/Header"].attrs["NumPart_Total"][4]

# Declare arrays for data
masses = zeros((n_parts,n_snapshots))
mass_from_AGB = zeros((n_parts,n_snapshots))
metallicity = zeros((n_parts,n_snapshots))
coord_parts = zeros((n_parts,3))
coord_sparts = zeros(3)
id = zeros((n_parts,n_snapshots))

id[:,0] = sim["/PartType0/ParticleIDs"]
coord_parts_unsorted = sim["/PartType0/Coordinates"]
coord_parts = sort_coords_by_id(id[:,0], coord_parts_unsorted)

# Read position of spart, if multiple sparts used ignore their positions
if n_sparts == 1:
	coord_sparts[:] = sim["/PartType4/Coordinates"][:]
	smoothing_length = sim["/PartType4/SmoothingLength"]
else:
	print("too many sparts for now")

# Read fields we are checking from snapshots and sort based on id
for i in range(n_snapshots):
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	print('reading snapshot '+str(i))
	id[:,i] = sim["/PartType0/ParticleIDs"]
	masses_unsorted = sim["/PartType0/Masses"]
	masses[:,i] = sort_by_id(id[:,i], masses_unsorted)
	mass_from_AGB_unsorted = sim["/PartType0/TotalMassFromAGB"]
	mass_from_AGB[:,i] = sort_by_id(id[:,i], mass_from_AGB_unsorted)
	metallicity_unsorted = sim["/PartType0/Metallicity"]
	metallicity[:,i] = sort_by_id(id[:,i], metallicity_unsorted)

# Now we can make some plots
figure()
subplot(111)

for i in range(100):
	#if ((distance(coord_parts[i,:],coord_sparts) < smoothing_length) & (metallicity[i,n_snapshots-1] != 0)):
	plot(metallicity[i,:],color='k',alpha=0.1, linewidth=0.5)

xlabel("snapshot")
ylabel("metallicity")
savefig("stellar_evolution_metallicity.png", dpi=200)

figure()
subplot(111)

for i in range(100):
	#if ((distance(coord_parts[i,:],coord_sparts) < smoothing_length) & (metallicity[i,n_snapshots-1] != 0)):
	plot(masses[i,:],color='k',alpha=0.1, linewidth=0.5)

xlabel("snapshot")
ylabel("mass")
savefig("stellar_evolution_mass.png", dpi=200)

figure()
subplot(111)

for i in range(n_parts):
	#if ((distance(coord_parts[i,:],coord_sparts) < smoothing_length) & (metallicity[i,n_snapshots-1] != 0)):
	plot(mass_from_AGB[i,:],color='k',alpha=0.1, linewidth=0.5)

xlabel("snapshot")
ylabel("mass from AGB")
savefig("stellar_evolution_mass_AGB.png", dpi=200)

total_part_mass = np.sum(masses,axis = 0)

figure()
subplot(111)

plot(total_part_mass/total_part_mass[0],color='k', linewidth=0.5)

xlabel("snapshot")
ylabel("total mass normalised by initial total mass")
savefig("stellar_evolution_total_mass.png", dpi=200)
