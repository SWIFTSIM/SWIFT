import numpy as np
import h5py as h5
from tqdm import tqdm
from numba import jit, prange

snapname = "eagle_0000/eagle_0000.hdf5"
fofname = "fof_output_0000/fof_output_0000.0.hdf5"
# snapname = "eagle_0000.hdf5"
# fofname = "fof_output_0000.hdf5"

######################################################

snap = h5.File(snapname, "r")

nogrp_grp_id = int(snap["/Parameters"].attrs.get("FOF:group_id_default"))

pos_gas = snap["/PartType0/Coordinates"][:, :]
ids_gas = snap["/PartType0/ParticleIDs"][:]
grp_gas = snap["/PartType0/FOFGroupIDs"][:]
mass_gas = snap["/PartType0/Masses"][:]

pos_DM = snap["/PartType1/Coordinates"][:, :]
ids_DM = snap["/PartType1/ParticleIDs"][:]
grp_DM = snap["/PartType1/FOFGroupIDs"][:]
mass_DM = snap["/PartType1/Masses"][:]

pos_star = snap["/PartType4/Coordinates"][:, :]
ids_star = snap["/PartType4/ParticleIDs"][:]
grp_star = snap["/PartType4/FOFGroupIDs"][:]
mass_star = snap["/PartType4/Masses"][:]

####################################################

fof = h5.File(fofname, "r")
num_files = fof["/Header/"].attrs["NumFilesPerSnapshot"][0]
num_groups = fof["/Header/"].attrs["NumGroups_Total"][0]
fof.close()

fof_grp = np.zeros(num_groups, dtype=np.int32)
fof_size = np.zeros(num_groups, dtype=np.int32)
fof_mass = np.zeros(num_groups)

# Read the distributed catalog
offset = 0
for i in range(num_files):

    my_filename = fofname[:-6]
    my_filename = my_filename + str(i) + ".hdf5"
    fof = h5.File(my_filename, "r")

    my_fof_grp = fof["/Groups/GroupIDs"][:]
    my_fof_size = fof["/Groups/Sizes"][:]
    my_fof_mass = fof["/Groups/Masses"][:]

    num_this_file = fof["/Header"].attrs["NumGroups_ThisFile"][0]
    fof.close()

    fof_grp[offset : offset + num_this_file] = my_fof_grp
    fof_size[offset : offset + num_this_file] = my_fof_size
    fof_mass[offset : offset + num_this_file] = my_fof_mass

    offset += num_this_file

####################################################

boxsize = snap["/Header"].attrs.get("BoxSize")[0]
N_DM = snap["/Header"].attrs.get("NumPart_ThisFile")[1]

l = 0.2 * boxsize / float(N_DM) ** (1.0 / 3.0)

print("Checking snapshot :", snapname)
print("Checking catalogue:", fofname)
print("L:", boxsize)
print("N_DM:", N_DM)
print("Linking length:", l)
print("")

####################################################


@jit(nopython=True, parallel=True, fastmath=True)
def nearest(dx, L=boxsize):
    mask1 = dx > 0.5 * L
    mask2 = dx < -0.5 * L
    if np.sum(mask1):
        dx[mask1] = dx[mask1] - L
    if np.sum(mask2):
        dx[mask2] = dx[mask2] + L
    return dx


####################################################

# Verify the content of the catalog
num_groups = np.size(fof_grp)
print("Catalog has", num_groups, "groups")


def check_fof_size(i):
    my_grp = fof_grp[i]
    my_size = fof_size[i]

    mask_gas = grp_gas == my_grp
    mask_DM = grp_DM == my_grp
    mask_star = grp_star == my_grp

    total = np.sum(mask_gas) + np.sum(mask_DM) + np.sum(mask_star)

    if total != my_size:
        print(
            "Grp",
            my_grp,
            "has size=",
            my_size,
            "but",
            total,
            "particles in the snapshot",
        )
        exit()


for i in tqdm(range(num_groups)):
    check_fof_size(i)

print("All group sizes match the particles")
####################################################

# Verify group masses
num_groups = np.size(fof_grp)
print("Catalog has", num_groups, "groups")


def check_fof_masses(i):
    my_grp = fof_grp[i]
    my_mass = fof_mass[i]

    mask_gas = grp_gas == my_grp
    mask_DM = grp_DM == my_grp
    mask_star = grp_star == my_grp

    total = (
        np.sum(mass_gas[mask_gas])
        + np.sum(mass_DM[mask_DM])
        + np.sum(mass_star[mask_star])
    )

    ratio = total / my_mass

    if ratio > 1.0001 or ratio < 0.9999:
        print(
            "Grp",
            my_grp,
            "has mass=",
            my_mass,
            "but particles in the snapshot have mass",
            total,
        )
        exit()


for i in tqdm(range(num_groups)):
    check_fof_masses(i)

print("All group masses match the particles")
####################################################

# Test the stand-alone stars
mask = grp_star == nogrp_grp_id
num_stars = np.sum(mask)
print("Found %d stars not in groups" % num_stars)
my_pos_star = pos_star[mask, :]
my_ids_star = ids_star[mask]
my_grp_star = grp_star[mask]
my_pos_DM = pos_DM[:, :]
my_ids_DM = ids_DM[:]
my_grp_DM = grp_DM[:]

# @jit(nopython=True, parallel=True, fastmath=True)
def check_stand_alone_star(i):
    pos = my_pos_star[i, :]
    grp = my_grp_star[i]

    dx = pos[0] - my_pos_DM[:, 0]
    dy = pos[1] - my_pos_DM[:, 1]
    dz = pos[2] - my_pos_DM[:, 2]

    dx = nearest(dx)
    dy = nearest(dy)
    dz = nearest(dz)

    # Identify the nearest DM particle
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    select = np.argmin(r2)

    # If the nearest DM particle is in a group --> mistake
    target_grp = my_grp_DM[select]
    if target_grp != nogrp_grp_id and r2[select] < l * l:
        print("Found a star without group whose nearest DM particle is in a group!")
        print("Star: id=", my_ids_star[i], "pos=", pos, "grp=", grp)
        print(
            "DM: id=",
            my_ids_DM[select],
            "pos=",
            my_pos_DM[select],
            "grp=",
            my_grp_DM[select],
        )
        print("r=", np.sqrt(r2[select]))
        # exit()


for i in tqdm(range(num_stars)):
    check_stand_alone_star(i)

print("All stand-alone stars OK!")

####################################################

# Test the stars in groups
mask = grp_star != nogrp_grp_id
num_stars = np.sum(mask)
print("Found %d stars in groups" % num_stars)
my_pos_star = pos_star[mask, :]
my_ids_star = ids_star[mask]
my_grp_star = grp_star[mask]
my_pos_DM = pos_DM[:, :]
my_ids_DM = ids_DM[:]
my_grp_DM = grp_DM[:]


@jit(nopython=True, parallel=True, fastmath=True)
def test_stars_in_group(i):
    pos = my_pos_star[i, :]
    grp = my_grp_star[i]

    dx = pos[0] - my_pos_DM[:, 0]
    dy = pos[1] - my_pos_DM[:, 1]
    dz = pos[2] - my_pos_DM[:, 2]

    dx = nearest(dx)
    dy = nearest(dy)
    dz = nearest(dz)

    # Identify the nearest DM particle
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    select = np.argmin(r2)

    # If the nearest DM particle is not in the same group --> mistake
    target_grp = my_grp_DM[select]
    if target_grp != grp and r2[select] < l * l:
        print(
            "Found a star in a group whose nearest DM particle is in a different group!"
        )
        print("Star: id=", my_ids_star[i], "pos=", pos, "grp=", grp)
        print(
            "DM: id=", my_ids_DM[select], "pos=", my_pos_DM[select], "grp=", target_grp
        )
        print("r=", np.sqrt(r2[select]))
        # exit()


for i in tqdm(range(num_stars)):
    test_stars_in_group(i)

print("All stars in groups OK!")

####################################################

# Test the stand-alone gas
mask = grp_gas == nogrp_grp_id
num_gas = np.sum(mask)
print("Found %d gas not in groups" % num_gas)
my_pos_gas = pos_gas[mask, :]
my_ids_gas = ids_gas[mask]
my_grp_gas = grp_gas[mask]
my_pos_DM = pos_DM[:, :]
my_ids_DM = ids_DM[:]
my_grp_DM = grp_DM[:]


@jit(nopython=True, parallel=True, fastmath=True)
def test_stand_alone_gas(i):
    pos = my_pos_gas[i, :]
    grp = my_grp_gas[i]

    dx = pos[0] - my_pos_DM[:, 0]
    dy = pos[1] - my_pos_DM[:, 1]
    dz = pos[2] - my_pos_DM[:, 2]

    dx = nearest(dx)
    dy = nearest(dy)
    dz = nearest(dz)

    # Identify the nearest DM particle
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    select = np.argmin(r2)

    # If the nearest DM particle is in a group --> mistake
    target_grp = my_grp_DM[select]
    if target_grp != nogrp_grp_id and r2[select] < l * l:
        print("Found a gas without group whose nearest DM particle is in a group!")
        print("Gas: id=", my_ids_gas[i], "pos=", pos, "grp=", grp)
        print(
            "DM: id=",
            my_ids_DM[select],
            "pos=",
            my_pos_DM[select],
            "grp=",
            my_grp_DM[select],
        )
        print("r=", np.sqrt(r2[select]))
        # exit()


for i in tqdm(range(num_gas)):
    test_stand_alone_gas(i)

print("All stand-alone gas OK!")

####################################################

# Test the gas in groups
mask = grp_gas != nogrp_grp_id
num_gas = np.sum(mask)
print("Found %d gas in groups" % num_gas)
my_pos_gas = pos_gas[mask, :]
my_ids_gas = ids_gas[mask]
my_grp_gas = grp_gas[mask]
my_pos_DM = pos_DM[:, :]
my_ids_DM = ids_DM[:]
my_grp_DM = grp_DM[:]


@jit(nopython=True, parallel=True, fastmath=True)
def test_gas_in_groups(i):
    pos = my_pos_gas[i, :]
    grp = my_grp_gas[i]

    dx = pos[0] - my_pos_DM[:, 0]
    dy = pos[1] - my_pos_DM[:, 1]
    dz = pos[2] - my_pos_DM[:, 2]

    dx = nearest(dx)
    dy = nearest(dy)
    dz = nearest(dz)

    # Identify the nearest DM particle
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    select = np.argmin(r2)

    # If the nearest DM particle is not in the same group --> mistake
    target_grp = my_grp_DM[select]
    if target_grp != grp and r2[select] < l * l:
        print(
            "Found a gas in a group whose nearest DM particle is in a different group!"
        )
        print("Gas: id=", my_ids_gas[i], "pos=", pos, "grp=", grp)
        print(
            "DM: id=", my_ids_DM[select], "pos=", my_pos_DM[select], "grp=", target_grp
        )
        print("r=", np.sqrt(r2[select]))
        # exit()


for i in tqdm(range(num_gas)):
    test_gas_in_groups(i)

print("All gas in groups OK!")

####################################################

# Test the stand-alone DM
mask = grp_DM == nogrp_grp_id
num_DM = np.sum(mask)
print("Found %d DM not in groups" % num_DM)
my_pos_DM = pos_DM[mask, :]
my_ids_DM = ids_DM[mask]
my_grp_DM = grp_DM[mask]


@jit(nopython=True, parallel=True, fastmath=True)
def test_stand_alone_DM(i):
    pos = my_pos_DM[i, :]
    grp = my_grp_DM[i]

    dx = pos[0] - pos_DM[:, 0]
    dy = pos[1] - pos_DM[:, 1]
    dz = pos[2] - pos_DM[:, 2]

    dx = nearest(dx)
    dy = nearest(dy)
    dz = nearest(dz)

    # Identify the nearest DM particle
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    mask = np.logical_and(r2 < l * l, r2 > 0.0)

    # If the nearest DM particle is in a group --> mistake
    if not np.all(grp_DM[mask] == nogrp_grp_id):
        print("Found a DM without group with some DM particle within l in a group!")
        print("DM:    id=", my_ids_DM[i], "pos=", pos, "grp=", grp)
        for j in range(np.sum(mask)):
            if grp_DM[mask][j] != nogrp_grp_id:
                print(
                    "Other: id=",
                    ids_DM[mask][j],
                    "pos=",
                    pos_DM[mask, :][j, :],
                    "grp=",
                    grp_DM[mask][j],
                    "r=",
                    np.sqrt(r2[mask][j]),
                )


for i in tqdm(range(num_DM)):
    test_stand_alone_DM(i)

print("All stand-alone DM OK!")

####################################################

# Test the DM in groups
mask = grp_DM != nogrp_grp_id
num_DM = np.sum(mask)
print("Found %d DM in groups" % num_DM)
my_pos_DM = pos_DM[mask, :]
my_ids_DM = ids_DM[mask]
my_grp_DM = grp_DM[mask]


@jit(nopython=True, parallel=True, fastmath=True)
def test_DM_in_groups(i):
    pos = my_pos_DM[i, :]
    grp = my_grp_DM[i]

    dx = pos[0] - pos_DM[:, 0]
    dy = pos[1] - pos_DM[:, 1]
    dz = pos[2] - pos_DM[:, 2]

    dx = nearest(dx)
    dy = nearest(dy)
    dz = nearest(dz)

    # Identify the nearest DM particle
    r2 = dx ** 2 + dy ** 2 + dz ** 2
    mask = r2 < l * l

    # If the nearest DM particle is not in the same group --> mistake
    if not np.all(grp_DM[mask] == grp):
        print(
            "Found a DM in a group whose DM particles within l are in a different group!"
        )
        print("DM:    id=", my_ids_DM[i], "pos=", pos, "grp=", grp)
        for j in range(np.sum(mask)):
            if grp_DM[mask][j] != grp:
                print(
                    "Other: id=",
                    ids_DM[mask][j],
                    "pos=",
                    pos_DM[mask, :][j, :],
                    "grp=",
                    grp_DM[mask][j],
                    "r=",
                    np.sqrt(r2[mask][j]),
                )


for i in tqdm(range(num_DM)):
    test_DM_in_groups(i)

print("All DM in groups OK!")
