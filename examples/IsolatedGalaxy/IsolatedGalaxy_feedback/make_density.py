import h5py
import sys
import numpy as np
import os

# Load density from ic folder (e.g. a SPHENIX run with the same setup)
initial_density_folder = sys.argv[2]
snapshot_path_initial_density = os.path.join(initial_density_folder, "output_0000.hdf5")

def write_or_create(group, name, data, dtype):
    if name in group:
        group[name][...] = data
    else:
        group.create_dataset(name, data=data, dtype=dtype)

with h5py.File(snapshot_path_initial_density, "r") as f:
    ids_ref = f["/PartType0/ParticleIDs"][:]
    rho_init = f["/PartType0/Densities"][:] # add [:] to fully load the data into memory, otherwise it is closed after the file is closed
    masses = f["/PartType0/Masses"][:]
    InternalEnergy = f["/PartType0/InternalEnergies"][:]

    # For python indexing:
    ids_ref = ids_ref - 1

# Create path to lowres8.hdf5 file (or other init files)
output_dir = sys.argv[1]
ic_path = os.path.join(output_dir, "fid.hdf5")

with h5py.File(ic_path, "r+") as f:
    N = f["PartType0/ParticleIDs"].shape[0]

    # Material IDs
    mat_ids = np.zeros(N)

    # Sort the densities of SPHENIX to the order of the reference (0 to N) -> that way the correct densities are assigned
    order = np.argsort(ids_ref)    
    rho_init = rho_init[order]

    # Also do this for the masses and internal energies (although these are roughly equal for all particles)
    masses = masses[order]
    InternalEnergy = InternalEnergy[order]

    # Particle group
    grp = f["PartType0"]
    write_or_create(grp, "Masses", masses, "f")
    write_or_create(grp, "InternalEnergies", InternalEnergy, "f")
    write_or_create(grp, "Density", rho_init, "f")
    write_or_create(grp, "MaterialIDs", mat_ids, "L")

print("Datasets added successfully.")