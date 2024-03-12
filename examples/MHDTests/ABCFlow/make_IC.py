#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters
rho = 1.0
cs2 = 3025.0
gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))
Bi_fraction = 1e-3

# output file
fileOutputName = "ABCFlow.hdf5"

def open_IAfile(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLength"][:]
    return pos, h

def add_other_particle_properties(pos,h,a,b,c,V0,kb,kv,Vz_factor,L):
    vol = L ** 3 
    N = len(h)

    # initializing arrays with particle properties
    v = np.zeros((N, 3))
    B = np.zeros((N, 3))
    A = np.zeros((N, 3))
    ids = np.linspace(1, N, N)
    m = np.ones(N) * rho * vol / N
    u = np.ones(N) * u0
    
    # setting constants
    kv0 = 2 * np.pi / L * kv
    kb0 = 2 * np.pi / L * kb
    Beq0 = np.sqrt(rho) * V0
    B0 = Bi_fraction * Beq0
    Norm = 1/np.sqrt(a**2+b**2+c**2)

    # rescaling the box to size L
    pos *= L

    # setting the velocity profile
    v[:,0] = (a * np.sin(kv0 * pos[:,2])+c * np.cos(kv0 * pos[:,1]))
    v[:,1] = (b * np.sin(kv0 * pos[:,0])+a * np.cos(kv0 * pos[:,2]))
    v[:,2] = (c * np.sin(kv0 * pos[:,1])+b * np.cos(kv0 * pos[:,0]))
    v[:,2] *= Vz_factor
    v*=V0 * Norm

    # setting the initial magnetic field
    # main mode for A.B. formula 6 from 1206.5186

    B[:, 0] = -(np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1]))
    B[:, 1] = -(np.cos(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2]))
    B[:, 2] = -(np.cos(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0]))
    B *= B0

    A[:, 0] = np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1])
    A[:, 1] = np.sin(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2])
    A[:, 2] = np.sin(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0])
    A0 = B0 / kb0
    A *= A0
    
    return v,B,A,ids,m,u


if __name__ == "__main__":

    import argparse as ap

    parser = ap.ArgumentParser(description="Generate ICs for ABC flow")

    parser.add_argument(
        "-V",
        "--rms_velocity",
        help="root mean square velocity of the flow",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "-P",
        "--IA_path",
        help="path to particle itinial arrangement file",
        default='./IAfiles/glassCube_32.hdf5',
        type=str,
    )
    parser.add_argument(
        "-kv",
        "--velocity_wavevector",
        help="wavelength of the velocity field",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-kb",
        "--magnetic_wavevector",
        help="wavelength of the initial magnetic field",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-L",
        "--boxsize",
        help="dimensions of the simulation box",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "-z",
        "--Vz_factor",
        help="multiplyier for velocity in z direciton",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "-A",
        "--A",
        help="constant A of the ABC flow",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "-B",
        "--B",
        help="constant B of the ABC flow",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "-C",
        "--C",
        help="constant C of the ABC flow",
        default=1.0,
        type=float,
    )

    args = parser.parse_args()
    pos,h = open_IAfile(args.IA_path)
    v,B,A,ids,m,u = add_other_particle_properties(pos,h,args.A,args.B,args.C,args.rms_velocity,args.magnetic_wavevector,args.velocity_wavevector,args.Vz_factor,args.boxsize)
    L=args.boxsize
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
        print('... done')
    except Exception as e: print('Error: '+str(e))
