################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2021 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
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
################################################################################

import h5py
import numpy as np

# Generates a swift IC file for the OrzagTang vortex in a periodic box


#fileOutputName = "stacked_28.hdf5"

# ---------------------------------------------------

Npside = 8
glass = h5py.File("glassCube_"+str(Npside)+".hdf5", "r")
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

Nold = len(h)

print(Nold)

times = 4
fileOutputName = "stacked_"+str(Npside*times)+".hdf5"

############ NEW
cx = times
cy = times
cz = times

lx = 1.0
ly = 1.0
lz = 1.0 

pnew = np.zeros((int(Nold * cx * cy * cz), 3))
hnew = np.zeros(int(Nold * cx * cy * cz))
N = Nold * cx * cy * cz

k = 0
c0 = 0
c1 = Nold
for i in range(0, cx):
    for j in range(0, cy):
        for k in range(0, cz):
            print(i, j, k)
            pnew[c0:c1, 0] = pos[:, 0] + i
            pnew[c0:c1, 1] = pos[:, 1] + j
            pnew[c0:c1, 2] = pos[:, 2] + k
            hnew[c0:c1] = h[:]
            c0 = c0 + Nold
            c1 = c1 + Nold

#print(len(pnew[:, 1]), " / ", N)
pos = pnew
pnew = 0
pos[:, 0] = pos[:, 0] * lx / cx  # -0.5
pos[:, 1] = pos[:, 1] * ly / cy  # -0.5
pos[:, 2] = pos[:, 2] * lz / cz
h = hnew / cx
hnew = 0
vol = 1.0
vol = lx * ly * lz

print(len(h))

fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [lx, ly, lz]  #####
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("SmoothingLength", data=h, dtype="f")

fileOutput.close()
