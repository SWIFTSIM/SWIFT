###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
##############################################################################

import h5py
import numpy as np

# Parameters
fileOutputName = "DecayingTurbulence.hdf5"

# Parameters
rho0 = 1.0
P = 1.0
cs2 = 3025.0
gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))
vrms = 1
Brms = 1


# smoothing kernel optimized for numpy Price 1012.1885 (6)
def open_IAfile(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLength"][:]
    return pos, h


def add_other_particle_properties(
    pos, h, L
):
    # size and number of particles
    vol = L ** 3
    N = len(h)

    # initializing arrays with particle properties
    v = np.zeros((N, 3))
    B = np.zeros((N, 3))
    A = np.zeros((N, 3))
    ids = np.linspace(1, N, N)
    m = np.ones(N) * rho0 * vol / N
    u = np.ones(N) * u0

    # rescaling the box to size L
    pos *= L
    h *= L

    # seeding turbulence (helical)
    v0 = vrms
    kv0 = 2*np.pi/L*10

    v[:, 0] += -(np.sin(kv0 * pos[:, 2]) - np.cos(kv0 * pos[:, 1]))
    v[:, 1] += -(np.cos(kv0 * pos[:, 0]) - np.cos(kv0 * pos[:, 2]))
    v[:, 2] += -(np.cos(kv0 * pos[:, 1]) - np.cos(kv0 * pos[:, 0]))
    v *= v0

    # seeding magnetic field (helical)
    B0 = Brms
    kb0 = 2*np.pi/L*1

    B[:, 0] += -(np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1]))
    B[:, 1] += -(np.cos(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2]))
    B[:, 2] += -(np.cos(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0]))
    B *= B0

    A[:, 0] += np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1])
    A[:, 1] += np.sin(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2])
    A[:, 2] += np.sin(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0])
    A0 = B0 / kb0
    A *= A0
 
    return pos, h, v, B, A, ids, m, u


if __name__ == "__main__":

    import argparse as ap

    parser = ap.ArgumentParser(description="Generate ICs for decaying turbulence box")

    parser.add_argument(
        "-L",
        "--boxsize",
        help="dimensions of the simulation box",
        default=2 * np.pi,
        type=float,
    )
    parser.add_argument(
        "-P",
        "--IA_path",
        help="path to particle itinial arrangement file",
        default="./IAfiles/glassCube_32.hdf5",
        type=str,
    )
 
    args = parser.parse_args()
    pos, h = open_IAfile(args.IA_path)
    pos, h, v, B, A, ids, m, u = add_other_particle_properties(
        pos,
        h,
        args.boxsize,
    )
    L = args.boxsize
    N = len(h)
    # File
    try:
        fileOutput = h5py.File(fileOutputName, "w")
        print("Wrinting ICs ...")
        # Header
        grp = fileOutput.create_group("/Header")
        grp.attrs["BoxSize"] = [L, L, L]  #####
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
        grp.create_dataset("Velocities", data=v, dtype="f")
        grp.create_dataset("Masses", data=m, dtype="f")
        grp.create_dataset("SmoothingLength", data=h, dtype="f")
        grp.create_dataset("InternalEnergy", data=u, dtype="f")
        grp.create_dataset("ParticleIDs", data=ids, dtype="L")
        grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
        grp.create_dataset("MagneticVectorPotentials", data=A, dtype="f")
        fileOutput.close()
        print("... done")
    except Exception as e:
        print("Error: " + str(e))
