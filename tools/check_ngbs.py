#!/usr/bin/env python

import h5py as h
import numpy as np
import matplotlib

matplotlib.use("Agg")
from pylab import *
import os.path

kernel_gamma = 1.825742
kernel_gamma2 = kernel_gamma * kernel_gamma
kernel_gamma_dim = np.power(kernel_gamma, 3)
hydro_dimension_unit_sphere = 4.0 * np.pi / 3.0
kernel_norm = hydro_dimension_unit_sphere * kernel_gamma_dim
error = False

inputFile1 = ""
inputFile2 = ""

# Compare the values of two floats
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


# Check list of density neighbours and check that they are correct.
def check_density_neighbours(
    pids, ngb_ids_naive, ngb_ids_sort, mask, pos, h_naive, h_sort, num_invalid, acc
):

    for k in range(0, num_invalid):

        # Filter neighbour lists for valid particle ids
        filter_neigh_naive = [i for i in ngb_ids_naive[mask][k] if i > -1]
        filter_neigh_sort = [i for i in ngb_ids_sort[mask][k] if i > -1]

        # Check neighbour lists for differences
        id_list = set(filter_neigh_naive).symmetric_difference(set(filter_neigh_sort))

        # Check for duplicate IDs
        duplicate_check_naive = len(filter_neigh_naive) != len(set(filter_neigh_naive))
        duplicate_check_sort = len(filter_neigh_sort) != len(set(filter_neigh_sort))

        if duplicate_check_naive:
            print("Duplicate neighbour ID found in: ", inputFile1)
            print(filter_neigh_naive)
            return True

        if duplicate_check_sort:
            print("Duplicate neighbour ID found in: ", inputFile2)
            print(filter_neigh_sort)
            return True

        pid = pids[mask][k]

        # Loop over discrepancies and check if they are actually neighbours
        for pjd in id_list:
            pi_pos = pos[np.where(pids == pid)]
            pj_pos = pos[np.where(pids == pjd)]

            hi = h_naive[np.where(pids == pid)]

            dx = pi_pos[0][0] - pj_pos[0][0]
            dy = pi_pos[0][1] - pj_pos[0][1]
            dz = pi_pos[0][2] - pj_pos[0][2]

            # Correct for BCs
            dx = nearest(dx)
            dy = nearest(dy)
            dz = nearest(dz)

            r2 = dx * dx + dy * dy + dz * dz

            hig2 = hi * hi * kernel_gamma2

            diff = abs(r2 - hig2)

            print(
                "Particle {} is missing {}, hig2: {}, r2: {}, |r2 - hig2|: {}".format(
                    pid, pjd, hig2, r2, diff
                )
            )

            if diff < acc * hig2:
                print("Missing interaction due to precision issue will be ignored.")
            else:
                hi_2 = h_sort[np.where(pids == pid)]

                # If a neigbour is missing and the particle has the same h throw
                # an error.
                if isclose(hi, hi_2):
                    print(
                        "Missing interaction found but particle has the same smoothing length (hi_1: %e, hi_2: %e)."
                        % (hi, hi_2)
                    )
                    return True
                else:
                    print(
                        "Missing interaction due to different smoothing lengths will be ignored (hi_1: %e, hi_2: %e)."
                        % (hi, hi_2)
                    )

    return False


# Check list of force neighbours and check that they are correct.
def check_force_neighbours(
    pids, ngb_ids_naive, ngb_ids_sort, mask, pos, h_naive, h_sort, num_invalid, acc
):

    error_val = False

    for k in range(0, num_invalid):

        # Filter neighbour lists for valid particle ids
        filter_neigh_naive = [i for i in ngb_ids_naive[mask][k] if i > -1]
        filter_neigh_sort = [i for i in ngb_ids_sort[mask][k] if i > -1]

        # Check neighbour lists for differences
        id_list = set(filter_neigh_naive).symmetric_difference(set(filter_neigh_sort))

        pid = pids[mask][k]

        # Loop over discrepancies and check if they are actually neighbours
        for pjd in id_list:
            pi_pos = pos[np.where(pids == pid)]
            pj_pos = pos[np.where(pids == pjd)]

            hi = h_naive[np.where(pids == pid)]
            hj = h_naive[np.where(pids == pjd)]

            dx = pi_pos[0][0] - pj_pos[0][0]
            dy = pi_pos[0][1] - pj_pos[0][1]
            dz = pi_pos[0][2] - pj_pos[0][2]

            # Correct for BCs
            dx = nearest(dx)
            dy = nearest(dy)
            dz = nearest(dz)

            r2 = dx * dx + dy * dy + dz * dz

            hig2 = hi * hi * kernel_gamma2
            hjg2 = hj * hj * kernel_gamma2

            diff = abs(r2 - max(hig2, hjg2))

            print(
                "Particle {} is missing {}, hig2: {}, hjg2: {}, r2: {}, |r2 - max(hig2,hjg2)|: {}".format(
                    pid, pjd, hig2, hjg2, r2, diff
                )
            )

            if diff < acc * max(hig2, hjg2):
                print("Missing interaction due to precision issue will be ignored.")
            else:
                hi_2 = h_sort[np.where(pids == pid)]
                if isclose(hi, hi_2):
                    print(
                        "Missing interaction due to the same smoothing lengths will not be ignored (hi_1: %e, hi_2: %e)."
                        % (hi, hi_2)
                    )
                    error_val = True
                else:
                    print(
                        "Missing interaction due to different smoothing lengths will be ignored (hi_1: %e, hi_2: %e)."
                        % (hi, hi_2)
                    )

    return error_val


def nearest(dx):
    if dx > 0.5 * box_size:
        return dx - box_size
    elif dx < -0.5 * box_size:
        return dx + box_size
    else:
        return dx


# Parse command line arguments
if len(sys.argv) < 3:
    print("Error: pass input files as arguments")
    sys.exit()
else:
    inputFile1 = sys.argv[1]
    inputFile2 = sys.argv[2]
    if os.path.exists(inputFile1) != 1:
        print("\n{} does not exist!\n".format(inputFile1))
        sys.exit()
    if os.path.exists(inputFile2) != 1:
        print("\n{} does not exist!\n".format(inputFile2))
        sys.exit()

# Open input files
file_naive = h.File(inputFile1, "r")
file_sort = h.File(inputFile2, "r")

box_size = file_naive["/Header"].attrs["BoxSize"][0]

# Read input file fields
ids_naive = file_naive["/PartType0/ParticleIDs"][:]
ids_sort = file_sort["/PartType0/ParticleIDs"][:]

h_naive = file_naive["/PartType0/SmoothingLength"][:]
h_sort = file_sort["/PartType0/SmoothingLength"][:]

pos_naive = file_naive["/PartType0/Coordinates"][:, :]
# pos_sort = file_sort["/PartType0/Coordinates"][:,:]

num_density_naive = file_naive["/PartType0/Num_ngb_density"][:]
num_density_sort = file_sort["/PartType0/Num_ngb_density"][:]

num_force_naive = file_naive["/PartType0/Num_ngb_force"][:]
num_force_sort = file_sort["/PartType0/Num_ngb_force"][:]

neighbour_ids_density_naive = file_naive["/PartType0/Ids_ngb_density"][:]
neighbour_ids_density_sort = file_sort["/PartType0/Ids_ngb_density"][:]

neighbour_ids_force_naive = file_naive["/PartType0/Ids_ngb_force"][:]
neighbour_ids_force_sort = file_sort["/PartType0/Ids_ngb_force"][:]


# wcount_naive = file_naive["/PartType0/Wcount"][:]
# wcount_sort = file_sort["/PartType0/Wcount"][:]
#
# wcount_naive = wcount_naive * np.power(h_naive,3) * kernel_norm
# wcount_sort = wcount_sort * np.power(h_sort,3) * kernel_norm

# Cross check
max_density_ngbs_naive = np.max(num_density_naive)
max_density_ngbs_sort = np.max(num_density_sort)
max_force_ngbs_naive = np.max(num_force_naive)
max_force_ngbs_sort = np.max(num_force_sort)

print("                   Min     Mean     Max ")
print("                   ---------------------")
print(
    "Ngbs density naiv: ",
    np.min(num_density_naive),
    np.mean(num_density_naive),
    max_density_ngbs_naive,
)
print(
    "Ngbs density sort: ",
    np.min(num_density_sort),
    np.mean(num_density_sort),
    max_density_ngbs_sort,
)
print(
    "Ngbs force naiv:   ",
    np.min(num_force_naive),
    np.mean(num_force_naive),
    max_force_ngbs_naive,
)
print(
    "Ngbs force sort:   ",
    np.min(num_force_sort),
    np.mean(num_force_sort),
    max_force_ngbs_sort,
)
# print "Wcount naiv:   ", np.min(wcount_naive), np.mean(wcount_naive), np.max(wcount_naive)
# print "Wcount sort:   ", np.min(wcount_sort), np.mean(wcount_sort), np.max(wcount_sort)

# Sort
index_naive = np.argsort(ids_naive)
index_sort = np.argsort(ids_sort)

num_density_naive = num_density_naive[index_naive]
num_density_sort = num_density_sort[index_sort]
num_force_naive = num_force_naive[index_naive]
num_force_sort = num_force_sort[index_sort]
ids_naive = ids_naive[index_naive]
ids_sort = ids_sort[index_sort]
neighbour_ids_density_naive = neighbour_ids_density_naive[index_naive]
neighbour_ids_density_sort = neighbour_ids_density_sort[index_sort]
neighbour_ids_force_naive = neighbour_ids_force_naive[index_naive]
neighbour_ids_force_sort = neighbour_ids_force_sort[index_sort]
# wcount_naive = wcount_naive[index_naive]
# wcount_sort = wcount_sort[index_sort]
h_naive = h_naive[index_naive]
h_sort = h_sort[index_sort]
pos_naive = pos_naive[index_naive]
# pos_sort = pos_sort[index_sort]

neighbour_length_naive = len(neighbour_ids_density_naive[0])
neighbour_length_sort = len(neighbour_ids_density_sort[0])

# Check that input files are logging the same number of neighbours
if neighbour_length_naive != neighbour_length_sort:
    print("Input files have logged different numbers of neighbour lengths!")
    print("{} has logged: {} neighbours".format(inputFile1, neighbour_length_naive))
    print("{} has logged: {} neighbours".format(inputFile2, neighbour_length_sort))
    exit(1)

if (
    max_density_ngbs_naive > neighbour_length_naive
    or max_force_ngbs_naive > neighbour_length_naive
    or max_density_ngbs_sort > neighbour_length_sort
    or max_force_ngbs_sort > neighbour_length_sort
):
    print("The number of neighbours has exceeded the number of neighbours logged.")
    print("Modify NUM_OF_NEIGHBOURS in hydro_part.h to log more neighbours.")
    print(
        "The highest neighbour count is: ",
        max(
            max_density_ngbs_naive,
            max_force_ngbs_naive,
            max_density_ngbs_sort,
            max_force_ngbs_sort,
        ),
    )
    exit(1)

# First check
print("\n                         Min    Max")
print("                         ----------")
print(
    "Differences for density:  ",
    min(num_density_naive - num_density_sort),
    max(num_density_naive - num_density_sort),
)
print(
    "Differences for force:    ",
    min(num_force_naive - num_force_sort),
    max(num_force_naive - num_force_sort),
)

# Get the IDs that are different
mask_density = num_density_naive != num_density_sort
mask_force = num_force_naive != num_force_sort
num_invalid_density = np.sum(mask_density)
num_invalid_force = np.sum(mask_force)

print("\nNum non-zero density: ", num_invalid_density)
print("Num non-zero force:   ", num_invalid_force)

print("\nParticle IDs with incorrect densities")
print("----------------------------------------")
print(ids_naive[mask_density])

# Check density neighbour lists
error += check_density_neighbours(
    ids_naive,
    neighbour_ids_density_naive,
    neighbour_ids_density_sort,
    mask_density,
    pos_naive,
    h_naive,
    h_sort,
    num_invalid_density,
    2e-6,
)

print("Num of density interactions", inputFile1)
print(num_density_naive[mask_density])

print("Num of density interactions", inputFile2)
print(num_density_sort[mask_density])

print("\nParticle IDs with incorrect forces")
print("------------------------------------")
print(ids_naive[mask_force])

# Check force neighbour lists
error += check_force_neighbours(
    ids_naive,
    neighbour_ids_force_naive,
    neighbour_ids_force_sort,
    mask_force,
    pos_naive,
    h_naive,
    h_sort,
    num_invalid_force,
    2e-6,
)

print("Num of force interactions", inputFile1)
print(num_force_naive[mask_force])

# print "Smoothing lengths", inputFile1
# print h_naive[mask_force]

print("Num of force interactions", inputFile2)
print(num_force_sort[mask_force])

# print "Smoothing lengths", inputFile2
# print h_sort[mask_force]

# Statistics of h difference
h_relative = (h_naive - h_sort) / h_naive
print(
    "h statistics: {} {} (Min, 1st Percentile)".format(
        np.min(h_relative), np.percentile(h_relative, 1)
    )
)
print(
    "h statistics: {} {} (Mean, Median)".format(
        np.mean(h_relative), np.median(h_relative)
    )
)
print(
    "h statistics: {} {} (Max, 99th Percentile)".format(
        np.max(h_relative), np.percentile(h_relative, 99)
    )
)

if error:
    print("\n------------------")
    print("Differences found.")
    print("------------------")
    exit(1)
else:
    print("\n---------------------")
    print("No differences found.")
    print("---------------------")
    exit(0)
