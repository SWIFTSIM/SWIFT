###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
##############################################################################

# Check the output done with swift

import h5py

filename = "testSelectOutput_0000.hdf5"
log_filename = "select_output.log"

# Read the simulation data
sim = h5py.File(filename, "r")
part0 = sim["/PartType0"]

# check presence / absence fields
if "Velocities" in part0:
    raise Exception("`Velocities` present in HDF5 but should not be written")

if "Coordinates" not in part0:
    raise Exception("`Coordinates` not present in HDF5 but should be written")

if "Masses" not in part0:
    raise Exception("`Masses` not present in HDF5 but should be written")

if "Densities" not in part0:
    raise Exception("`Densities` not present in HDF5 but should be written")


# check error detection
with open(log_filename, "r") as f:
    data = f.read()

if "Default:Masses_Gas" not in data:
    raise Exception("Input error in `Default:Masses_Gas` not detected")

if "Default:Pot_Gas" not in data:
    raise Exception("Parameter name error not detected for `SelectOutput:Pot_Gas`")
