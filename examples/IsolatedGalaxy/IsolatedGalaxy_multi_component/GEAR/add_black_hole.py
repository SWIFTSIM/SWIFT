#!/usr/bin/env python3
"""Seed a single central black hole in the multi-component galaxy ICs.

The downloaded galaxy_multi_component.hdf5 file has no PartType5 group. This
script adds one, placed at the gas centre of mass, so the example can
optionally be run with GEAR's BLACK_HOLES_GEAR module (see with_bh in
run.sh). It is idempotent: if the file already contains a PartType5 group,
it does nothing.

Usage:
    ./add_black_hole.py [galaxy_multi_component.hdf5]
"""

import sys

import h5py
import numpy as np

# Initial dynamical and subgrid mass of the seed black hole.
SEED_MASS_MSUN = 1.5e5
SOLAR_MASS_IN_G = 1.98841e33


def add_dataset(group, ref_group, name, value, name_out=None):
    """Create a length-1 dataset in `group` matching the dtype/shape of ref_group[name]."""
    ref = ref_group[name]
    data = np.asarray(value, dtype=ref.dtype)
    data = data.reshape((1,) + ref.shape[1:])
    group.create_dataset(name_out or name, data=data)


def main(filename):
    with h5py.File(filename, "r+") as f:
        if "PartType5" in f:
            print(f"{filename} already has a black hole particle, nothing to do.")
            return

        gas = f["PartType0"]
        gas_pos = gas["Coordinates"][:]
        gas_vel = gas["Velocities"][:]
        gas_mass = gas["Masses"][:]
        gas_h = gas["SmoothingLength"][:]

        com_pos = np.average(gas_pos, axis=0, weights=gas_mass)
        com_vel = np.average(gas_vel, axis=0, weights=gas_mass)
        # Initial guess only; the density loop iterates it to convergence.
        seed_h = np.mean(gas_h)

        max_id = 0
        for key in f:
            if key.startswith("PartType") and "ParticleIDs" in f[key]:
                max_id = max(max_id, int(np.max(f[key]["ParticleIDs"][:])))
        new_id = max_id + 1

        # Read from the file's own declared unit system (scalar in pNbody-written ICs, length-1 array in downloaded ones).
        unit_mass_in_g = np.ravel(f["Units"].attrs["Unit mass in cgs (U_M)"])[0]
        seed_mass = SEED_MASS_MSUN * SOLAR_MASS_IN_G / unit_mass_in_g

        bh = f.create_group("PartType5")
        add_dataset(bh, gas, "Coordinates", com_pos)
        add_dataset(bh, gas, "Velocities", com_vel)
        add_dataset(bh, gas, "Masses", [seed_mass])
        add_dataset(bh, gas, "ParticleIDs", [new_id])
        add_dataset(bh, gas, "SmoothingLength", [seed_h])
        # Subgrid mass starts equal to the dynamical mass, as there is no
        # separate resolved/unresolved distinction for a hand-placed seed.
        add_dataset(bh, bh, "Masses", [seed_mass], name_out="SubgridMasses")

        header = f["Header"]
        num_this_file = list(header.attrs["NumPart_ThisFile"])
        num_total = list(header.attrs["NumPart_Total"])
        num_this_file[5] = 1
        num_total[5] = 1
        header.attrs["NumPart_ThisFile"] = num_this_file
        header.attrs["NumPart_Total"] = num_total

    print(
        f"Added 1 black hole seed (id={new_id}, mass={SEED_MASS_MSUN:g} Msun) "
        f"at the gas centre of mass {com_pos} to {filename}."
    )


if __name__ == "__main__":
    filename = sys.argv[1] if len(sys.argv) > 1 else "galaxy_multi_component.hdf5"
    main(filename)
