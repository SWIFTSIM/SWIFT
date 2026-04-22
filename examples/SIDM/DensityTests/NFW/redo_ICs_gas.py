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

import h5py
import numpy as np
import sys


def cdm_to_SIDM(input_filename, output_filename, old_PT=1, new_PT=7):

    old_name = f"PartType{old_PT}"
    new_name = f"PartType{new_PT}"

    with h5py.File(input_filename, "r") as f_in, h5py.File(
        output_filename, "w"
    ) as f_out:

        # Copy all groups except the old DM group
        for key in f_in.keys():
            if key == old_name:
                f_in.copy(old_name, f_out, name=new_name)
            else:
                f_in.copy(key, f_out)

        # Copy and update header attributes
        header_in = f_in["Header"]
        header_out = f_out["Header"]

        ndm = len(f_in[old_name]["Coordinates"])

        for name in ["NumPart_ThisFile", "NumPart_Total"]:
            old_arr = np.array(header_in.attrs[name])
            new_arr = np.zeros(8, dtype=old_arr.dtype)
            new_arr[: len(old_arr)] = old_arr
            new_arr[old_PT] -= ndm
            new_arr[new_PT] += ndm
            header_out.attrs[name] = new_arr

        # Fix MassTable
        mt_old = np.array(header_in.attrs["MassTable"])
        mt = np.zeros(8, dtype=mt_old.dtype)
        mt[: len(mt_old)] = mt_old
        header_out.attrs["MassTable"] = mt

        # Fix NumPart_Total_HighWord
        mt_old = np.array(header_in.attrs["NumPart_Total_HighWord"])
        mt = np.zeros(8, dtype=mt_old.dtype)
        mt[: len(mt_old)] = mt_old
        header_out.attrs["NumPart_Total_HighWord"] = mt

        # the rest unchanged
        for attr in header_in.attrs:
            if attr not in [
                "NumPart_ThisFile",
                "NumPart_Total",
                "NumPart_Total_HighWord",
                "MassTable",
            ]:
                header_out.attrs[attr] = header_in.attrs[attr]

        # Add DM smoothing lengths
        rmax = 50
        h = 1.2348 * ((4 * np.pi * rmax ** 3) / 3 / ndm) ** (1.0 / 3.0)
        u = np.zeros(ndm)

        parttype7 = f_out["PartType0"]
        parttype7.create_dataset("SmoothingLength", data=np.repeat(h, ndm), dtype="f")
        parttype7.create_dataset("InternalEnergy", data=u, dtype="f")


def main():

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    cdm_to_SIDM(input_filename, output_filename, old_PT=1, new_PT=0)


if __name__ == "__main__":
    main()
