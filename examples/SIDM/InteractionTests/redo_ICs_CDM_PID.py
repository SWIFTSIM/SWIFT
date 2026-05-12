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


def fix_pids(input_filename, output_filename):
    with h5py.File(input_filename, "r") as f_in, h5py.File(
        output_filename, "w"
    ) as f_out:

        for key in f_in:
            if key != "PartType1":
                f_in.copy(key, f_out)
                continue

            grp_in = f_in["PartType1"]
            grp_out = f_out.create_group("PartType1")

            for name, ds in grp_in.items():
                if name == "ParticleIDs":
                    grp_out.create_dataset(name, data=ds[:] + 1)
                else:
                    f_in.copy(ds, grp_out)


def main():

    input_filename = "dSph_cusp.hdf5"
    output_filename = "dSph_cusp_fixed.hdf5"
    fix_pids(input_filename, output_filename)


if __name__ == "__main__":
    main()
