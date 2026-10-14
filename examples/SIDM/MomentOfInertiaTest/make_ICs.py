###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

# Converts the DM particles (PartType1) of EAGLE ICs into either gas
# particles (PartType0) or SIDM particles (PartType7). Every other particle
# type in the input is dropped, so the gas and SIDM ICs contain exactly the
# same particles (same positions, velocities, masses and IDs).

import h5py
import numpy as np
import sys
from scipy.spatial import cKDTree

DM_PT = 1
TARGET_PT = {"gas": 0, "sidm": 7}
GAS_INTERNAL_ENERGY = 1.0

# Neighbour number and kernel used for the initial smoothing length guess
# (48 neighbours with the cubic spline kernel, i.e. resolution_eta = 1.2348)
N_NGB = 48
KERNEL_GAMMA = 1.825742


def estimate_smoothing_lengths(pos, box_size, chunk=1_000_000):

    pos = np.mod(pos, box_size)
    tree = cKDTree(pos, boxsize=box_size)
    h = np.empty(len(pos), dtype=np.float32)
    for i in range(0, len(pos), chunk):
        r, _ = tree.query(pos[i : i + chunk], k=[N_NGB], workers=-1)
        h[i : i + chunk] = r[:, 0] / KERNEL_GAMMA
    return h


def dm_to(input_filename, output_filename, target):

    new_PT = TARGET_PT[target]
    old_name = f"PartType{DM_PT}"
    new_name = f"PartType{new_PT}"

    with (
        h5py.File(input_filename, "r") as f_in,
        h5py.File(output_filename, "w") as f_out,
    ):

        # Copy all non-particle groups (Header, Units, Cosmology, ...)
        for key in f_in.keys():
            if not key.startswith("PartType"):
                f_in.copy(key, f_out)

        header_in = f_in["Header"]
        header_out = f_out["Header"]
        dm = f_in[old_name]
        npart = dm["Coordinates"].shape[0]
        box_size = np.atleast_1d(header_in.attrs["BoxSize"])[0]

        # Copy the DM particle data to the new particle type, one dataset at a
        # time to avoid loading everything into memory at once
        grp = f_out.create_group(new_name)
        for field in ["Coordinates", "Velocities", "ParticleIDs"]:
            f_in.copy(dm[field], grp, name=field)

        if "Masses" in dm:
            f_in.copy(dm["Masses"], grp, name="Masses")
        else:
            mass = np.array(header_in.attrs["MassTable"])[DM_PT]
            grp.create_dataset("Masses", data=np.full(npart, mass, dtype=np.float32))

        # DM particles have no smoothing lengths: estimate them from the
        # distance to the N_NGB-th neighbour (identical for gas and SIDM).
        # A uniform guess is not enough, as SWIFT builds the top-level grid
        # from the initial h and particles in voids then outgrow their cells.
        h = estimate_smoothing_lengths(dm["Coordinates"][:], box_size)
        grp.create_dataset("SmoothingLength", data=h)

        if target == "gas":
            grp.create_dataset(
                "InternalEnergy",
                data=np.full(npart, GAS_INTERNAL_ENERGY, dtype=np.float32),
            )

        # Rewrite the particle counts: only the new particle type remains
        for name in ["NumPart_ThisFile", "NumPart_Total", "NumPart_Total_HighWord"]:
            old_arr = np.array(header_in.attrs[name])
            new_arr = np.zeros(8, dtype=old_arr.dtype)
            if name == "NumPart_Total_HighWord":
                new_arr[new_PT] = old_arr[DM_PT]
            else:
                new_arr[new_PT] = npart
            header_out.attrs[name] = new_arr

        # Masses are stored explicitly
        mt_old = np.array(header_in.attrs["MassTable"])
        header_out.attrs["MassTable"] = np.zeros(8, dtype=mt_old.dtype)

        if "NumFilesPerSnapshot" in header_in.attrs:
            header_out.attrs["NumFilesPerSnapshot"] = 1

    print(f"Wrote {npart} {target} particles to {output_filename}")


def main():

    if len(sys.argv) != 4 or sys.argv[3] not in TARGET_PT:
        sys.exit(1)

    dm_to(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
