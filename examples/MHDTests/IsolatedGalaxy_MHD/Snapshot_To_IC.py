###############################################################################


import h5py
import numpy as np
import matplotlib.pyplot as plt

# smoothing kernel optimized for numpy Price 1012.1885 (6)
def open_IAfile_gas(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    boxsize = IAfile["Header"].attrs["BoxSize"][0]
    t = IAfile["Header"].attrs["Time"][0]
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLengths"][:]
    v = IAfile["/PartType0/Velocities"][:]
    B = IAfile["/PartType0/MagneticFluxDensities"][:,:]
    m = IAfile["/PartType0/Masses"]
    u = IAfile["/PartType0/InternalEnergies"]
    pids = IAfile["/PartType0/ParticleIDs"]
    return boxsize,t,pos,h,v,B,m,u,pids

def open_IAfile_stars(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    boxsize = IAfile["Header"].attrs["BoxSize"][0]
    t = IAfile["Header"].attrs["Time"][0]
    pos = IAfile["/PartType4/Coordinates"][:, :]
    h = IAfile["/PartType4/SmoothingLengths"][:]
    v = IAfile["/PartType4/Velocities"][:]
    m = IAfile["/PartType4/Masses"]
    pids = IAfile["/PartType4/ParticleIDs"]
    return boxsize,t,pos,h,v,m,pids

# Make snapshot
def make_ICs_from_snapshot(path_to_snapshot,path_to_ICs="./CoolingHalo_new.hdf5"):
    ICfile = h5py.File(path_to_snapshot, "r")
    boxSize = ICfile["Header"].attrs["BoxSize"][0]
    Time = ICfile["Header"].attrs["Time"][0]
    print(Time)

    # Create the file
    filename = path_to_ICs
    file = h5py.File(filename, "w")
    grp = file.create_group("/Header")
    grp.attrs["BoxSize"] = boxSize
    #grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    #grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    #grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = Time
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["Dimension"] = 3
  
    # Read gas particle properties
    pos = ICfile["/PartType0/Coordinates"][:, :]
    h = ICfile["/PartType0/SmoothingLengths"][:]
    v = ICfile["/PartType0/Velocities"][:]
    B = ICfile["/PartType0/MagneticFluxDensities"][:,:]
    m = ICfile["/PartType0/Masses"]
    u = ICfile["/PartType0/InternalEnergies"]
    pids = ICfile["/PartType0/ParticleIDs"]
    Ng = len(pids)

    # Write gas particle properties
    grpG = file.create_group("/PartType0")

    ds = grpG.create_dataset("Coordinates", (Ng, 3), "d")
    ds[()] = pos
    ds = grpG.create_dataset("Velocities", (Ng, 3), "d")
    ds[()] = v
    ds = grpG.create_dataset("MagneticFluxDensities", (Ng, 3), "d")
    ds[()] = B
    ds = grpG.create_dataset("SmoothingLength", (Ng,), "f")
    ds[()] = h
    ds = grpG.create_dataset("Masses", (Ng,), "f")
    ds[()] = m
    ds = grpG.create_dataset("InternalEnergy", (Ng,), "f")
    ds[()] = u
    ds = grpG.create_dataset("ParticleIDs", (Ng,), "L")
    ds[()] = pids

    # Read star particle properties
    pos = ICfile["/PartType4/Coordinates"][:, :]
    h = ICfile["/PartType4/SmoothingLengths"][:]
    v = ICfile["/PartType4/Velocities"][:]
    m = ICfile["/PartType4/Masses"]
    pids = ICfile["/PartType4/ParticleIDs"]
    Ns = len(pids)

    # Create star particle properties
    grpS = file.create_group("/PartType4")
    ds = grpS.create_dataset("Coordinates", (Ns, 3), "d")
    ds[()] = pos
    ds = grpS.create_dataset("Velocities", (Ns, 3), "d")
    ds[()] = v
    ds = grpS.create_dataset("SmoothingLength", (Ns,), "f")
    ds[()] = h
    ds = grpS.create_dataset("Masses", (Ns,), "f")
    ds[()] = m
    ds = grpS.create_dataset("ParticleIDs", (Ns,), "L")
    ds[()] = pids

    Ntot = Ng+Ns
    grp.attrs["NumPart_Total"] = [Ng, 0, 0, 0, Ns, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [Ng, 0, 0, 0, Ns, 0]

    # Units
    grp = file.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 3.08567758E21
    grp.attrs["Unit mass in cgs (U_M)"] = 1.9891E43
    grp.attrs["Unit time in cgs (U_t)"] = 3.08567758E21/1E5
    grp.attrs["Unit current in cgs (U_I)"] = 2.088431e13
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    file.close()

make_ICs_from_snapshot("./EXP_CHASF05/CoolingHalo_0707.hdf5", "CoolingHaloNew.hdf5")
