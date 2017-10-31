import h5py as h
import numpy as np
import matplotlib
matplotlib.use("Agg")
from pylab import *
import os.path

kernel_gamma = 1.825742
kernel_gamma_dim = np.power(kernel_gamma,3)
hydro_dimension_unit_sphere = 4. * np.pi / 3.
kernel_norm = hydro_dimension_unit_sphere * kernel_gamma_dim

inputFile1 = ""
inputFile2 = ""

if len(sys.argv) < 3:
    print "Error: pass input files as arguments"
    sys.exit()
else:
    inputFile1 = sys.argv[1]
    inputFile2 = sys.argv[2]
    if os.path.exists(inputFile1) != 1:
        print "\n{} does not exist!\n".format(inputFile1)
        sys.exit()
    if os.path.exists(inputFile2) != 1:
        print "\n{} does not exist!\n".format(inputFile2)
        sys.exit()
    
file_naiv = h.File(inputFile1, "r")
file_sort = h.File(inputFile2, "r")

ids_naiv = file_naiv["/PartType0/ParticleIDs"][:]
ids_sort = file_sort["/PartType0/ParticleIDs"][:]

h_naiv = file_naiv["/PartType0/SmoothingLength"][:]
h_sort = file_sort["/PartType0/SmoothingLength"][:]

pos_naiv = file_naiv["/PartType0/Coordinates"][:,:]
pos_sort = file_sort["/PartType0/Coordinates"][:,:]

num_density_naiv = file_naiv["/PartType0/Num_ngb_density"][:]
num_density_sort = file_sort["/PartType0/Num_ngb_density"][:]

num_force_naiv = file_naiv["/PartType0/Num_ngb_force"][:]
num_force_sort = file_sort["/PartType0/Num_ngb_force"][:]

#wcount_naiv = file_naiv["/PartType0/Wcount"][:]
#wcount_sort = file_sort["/PartType0/Wcount"][:]
#
#wcount_naiv = wcount_naiv * np.power(h_naiv,3) * kernel_norm
#wcount_sort = wcount_sort * np.power(h_sort,3) * kernel_norm

# Cross check
print "Ngbs density naiv: ", np.min(num_density_naiv), np.mean(num_density_naiv), np.max(num_density_naiv)
print "Ngbs density sort: ", np.min(num_density_sort), np.mean(num_density_sort), np.max(num_density_sort)
print "Ngbs force naiv:   ", np.min(num_force_naiv), np.mean(num_force_naiv), np.max(num_force_naiv)
print "Ngbs force sort:   ", np.min(num_force_sort), np.mean(num_force_sort), np.max(num_force_sort)
#print "Wcount naiv:   ", np.min(wcount_naiv), np.mean(wcount_naiv), np.max(wcount_naiv)
#print "Wcount sort:   ", np.min(wcount_sort), np.mean(wcount_sort), np.max(wcount_sort)

# Sort
index_naiv = np.argsort(ids_naiv)
index_sort = np.argsort(ids_sort)

num_density_naiv = num_density_naiv[index_naiv]
num_density_sort = num_density_sort[index_sort]
num_force_naiv = num_force_naiv[index_naiv]
num_force_sort = num_force_sort[index_sort]
ids_naiv = ids_naiv[index_naiv]
ids_sort = ids_sort[index_sort]
#wcount_naiv = wcount_naiv[index_naiv]
#wcount_sort = wcount_sort[index_sort]
h_naiv = h_naiv[index_naiv]
h_sort = h_sort[index_sort]

# First check
print "Differences for density:  ", min(num_density_naiv - num_density_sort), max(num_density_naiv - num_density_sort)
print "Differences for force:    ", min(num_force_naiv - num_force_sort), max(num_force_naiv - num_force_sort)

# Get the IDs that are different
mask_density = num_density_naiv != num_density_sort
mask_force = num_force_naiv != num_force_sort

print "Num non-zero density: ", np.sum(mask_density)
print "Num non-zero force:   ", np.sum(mask_force)

print "\nParticle IDs with incorrect densities"
print "----------------------------------------"
print ids_naiv[mask_density]

print "Num of density interactions", inputFile1
print num_density_naiv[mask_density]

print "Num of density interactions", inputFile2
print num_density_sort[mask_density]

print "\nParticle IDs with incorrect forces"
print "------------------------------------"
print ids_naiv[mask_force]

print "Num of force interactions", inputFile1
print num_force_naiv[mask_force]

print "Num of force interactions", inputFile2
print num_force_sort[mask_force]

# Statistics of h difference
h_relative = (h_naiv - h_sort) / h_naiv
print "h statistics:", np.min(h_relative), np.percentile(h_relative, 1) 
print "h statistics:", np.mean(h_relative), np.median(h_relative) 
print "h statistics:", np.percentile(h_relative, 99) , np.max(h_relative) 

#print mean(h_naiv)

#exit()

# Make some plots of the wrong particles
#rcParams.update({'figure.figsize': (6,6)})
#figure()
#subplot(221)
#plot(pos_naiv[mask_force,0], pos_naiv[mask_force,1], 'bo', ms=2)
#for i in range(12):
#    plot([0, 1], [i / 12., i/12.], 'k-', lw=0.5, color='0.5')
#    plot([i / 12., i/12.], [0, 1], 'k-', lw=0.5, color='0.5')
#
#
#subplot(222)
#plot(pos_naiv[mask_force,0], pos_naiv[mask_force,2], 'bo', ms=2)
#for i in range(12):
#    plot([0, 1], [i / 12., i/12.], 'k-', lw=0.5, color='0.5')
#    plot([i / 12., i/12.], [0, 1], 'k-', lw=0.5, color='0.5')
#
#
#subplot(223)
#plot(pos_naiv[mask_force,1], pos_naiv[mask_force,2], 'bo', ms=2)
#for i in range(12):
#    plot([0, 1], [i / 12., i/12.], 'k-', lw=0.5, color='0.5')
#    plot([i / 12., i/12.], [0, 1], 'k-', lw=0.5, color='0.5')
#
#
#savefig("wrong_force.png")
