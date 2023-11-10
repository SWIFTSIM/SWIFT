#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters

r_in = 0.1
rho_0 = 1.0
P_in_0 = 10.0
P_out_0 = 0.1
B_0 = 1.0 
gamma = 5.0 / 3.0

fileOutputName = "MagneticBlastWave_LR.hdf5"

###---------------------------###

glass = h5py.File("glassCube_32.hdf5", "r")

unit_cell = glass["/PartType0/Coordinates"][:, :]
h_unit_cell = glass["/PartType0/SmoothingLength"][:]

N_unit_cell = len(h_unit_cell)
times = 2

###---------------------------###

cx = times
cy = times
cz = 1

lx = 2.0
ly = 2.0
lz = 2.0 / float(times)

lx_c = lx / 2
ly_c = ly / 2
lz_c = lz / 2

N = N_unit_cell * cx * cy * cz

pos = np.zeros((int(N), 3))
h = np.zeros(N)

c0 = 0
c1 = N_unit_cell

for i in range(0, cx):
    for j in range(0, cy):
        for k in range(0, cz):
            pos[c0:c1, 0] = unit_cell[:, 0] + i
            pos[c0:c1, 1] = unit_cell[:, 1] + j
            pos[c0:c1, 2] = unit_cell[:, 2] + k
            h[c0:c1] = h_unit_cell[:]
            c0 = c0 + N_unit_cell
            c1 = c1 + N_unit_cell

pos[:, 0] = pos[:, 0] * lx / cx
pos[:, 1] = pos[:, 1] * ly / cy
pos[:, 2] = pos[:, 2] * lz / cz
h = h * lx / cx

vol = lx * ly * lz

###---------------------------###

rot = np.sqrt((pos[:, 0] - lx_c) ** 2 + (pos[:, 1] - ly_c) ** 2)
#v = np.zeros((N, 3))
#B = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho_0 * vol / N
u = np.ones(N)
u[rot < r_in] = P_in_0 / (rho_0 * (gamma - 1))
u[rot >= r_in] = P_out_0 / (rho_0 * (gamma - 1))
#B[:, 0] = B_0


###---------------------------###
N2=int(2*N)
p=np.zeros((N2, 3))
p[:N,0]=pos[:,0]
p[N:,0]=pos[:,0]
p[:N,1]=pos[:,1]
p[N:,1]=pos[:,1]+1.0
p[:N,2]=pos[:,2]
p[N:,2]=pos[:,2]
pos=p
hh =np.zeros(N2)
hh[:N]=h
hh[N:]=h
h=hh
v=np.zeros((N2,3))
ids = np.linspace(1, N2, N2)
m = np.ones(N2) * rho_0 * vol / N_out
uu =np.zeros(N2)
uu[:N]=u
uu[N:]=u
u=uu
B = np.zeros((N2, 3))
A = np.zeros((N2, 3))
B[:N, 0] = B_0
B[N:, 0] = -B_0
A[:N, 2] = B_0 * pos[:N,1]
A[N:, 2] = B_0 * (2.0-pos[N:,1])

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [lx, 2.* ly, lz]  #####
grp.attrs["NumPart_Total"] = [2*N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [2*N, 0, 0, 0, 0, 0]
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
grp.create_dataset("MagneticVectorPotentials", data = A, dtype = 'f')

fileOutput.close()
