#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

# Parameters
rho = 1.0
cs2 = 3025.0
L = float(sys.argv[6])  # 1.0
kv = int(sys.argv[4])
kv0 = 2 * np.pi / L * kv
kb = int(sys.argv[5])
kb0 = 2 * np.pi / L * kb
V0 = float(sys.argv[1])  # 22.9
Vz_factor = float(sys.argv[7])
resistive_eta = float(sys.argv[2])
Flow_kind = int(sys.argv[8])
Beq0 = np.sqrt(rho) * V0
B0 = 1e-8 * Beq0
gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))

# output file
fileOutputName = "RobertsFlow.hdf5"

# If Make_New is true, the script makes new ICs from IA files and with specified velocity and magnetic field profile
# If false it uses selected snapshot to create IC file
Make_New = True
Continue = True

if Make_New:
    ###---------------------------###

    glass = h5py.File(sys.argv[3], "r")
    pos = glass["/PartType0/Coordinates"][:, :]
    h = glass["/PartType0/SmoothingLength"][:]

    N = len(h)
    vol = L ** 3

    ###---------------------------###

    v = np.zeros((N, 3))
    B = np.zeros((N, 3))
    A = np.zeros((N, 3))
    ids = np.linspace(1, N, N)
    m = np.ones(N) * rho * vol / N
    u = np.ones(N) * u0

    # rescaling the box to size L
    pos *= L

    # setting up flow

    #v[:, 0] = np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
    #v[:, 1] = -np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
    #v[:, 2] = Vz_factor * np.sqrt(2) * np.cos(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
   
    if Flow_kind==0:
        # setting up flow

        #v[:, 0] = np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
        #v[:, 1] = -np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
        #v[:, 2] = Vz_factor * np.sqrt(2) * np.cos(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])

        # velocity field consistent with A.B. runs 
        v[:, 0] = np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
        v[:, 1] = -np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
        v[:, 2] = Vz_factor * np.sqrt(2) * np.sin(kv0 * pos[:, 1]) * np.sin(kv0 * pos[:, 0])
        v *= V0

        # Exciting main RF mode

        #B[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        #B[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        #B[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        #B *= B0

        #A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        #A[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        #A[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        #A0 = B0 / kb0
        #A *= A0

        # main mode for A.B.

        B[:, 0] = -(np.sin(kb0 * pos[:, 2]) + np.sin(kb0 * pos[:, 1]))
        B[:, 1] = -(np.cos(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2]))
        B[:, 2] = -(np.cos(kb0 * pos[:, 1]) + np.sin(kb0 * pos[:, 0]))

        B *= B0

        A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.sin(kb0 * pos[:, 1])
        A[:, 1] = np.cos(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2])
        A[:, 2] = np.cos(kb0 * pos[:, 1]) + np.sin(kb0 * pos[:, 0])
        A0 = B0 / kb0
        A *= A0

        #random main mode for A.B.
        #B[:, 0] = 2*(np.random.rand(N)-0.5)
        #B[:, 1] = 2*(np.random.rand(N)-0.5)
        #B[:, 2] = 2*(np.random.rand(N)-0.5)
        #B *= B0

        #A[:, 0] = 2*(np.random.rand(N)-0.5)
        #A[:, 1] = 2*(np.random.rand(N)-0.5)
        #A[:, 2] = 2*(np.random.rand(N)-0.5)
        #A0 = B0 / kb0
        #A *= A0
    elif Flow_kind==1:
        # setting up flow (warning, this field has Vrms=sqrt(2))
        v[:, 0] = np.sin(kv0 * pos[:, 0])
        v[:, 1] = np.sin(kv0 * pos[:, 1])
        v[:, 2] = Vz_factor * (np.cos(kv0 * pos[:, 0]) - np.cos(kv0 * pos[:, 1]))
        v *= V0

        B[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        B[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        B[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        B *= B0

        A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        A[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        A[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        A0 = B0 / kb0
        A *= A0

    elif Flow_kind==2:
        # setting up flow (warning, this field has Vrms=sqrt(2))
        v[:, 0] = np.sin(kv0 * pos[:, 0])
        v[:, 1] = np.sin(kv0 * pos[:, 1])
        v[:, 2] = Vz_factor * (np.cos(kv0 * pos[:, 0]) + np.cos(kv0 * pos[:, 1]))
        v *= V0

        B[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        B[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        B[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        B *= B0

        A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        A[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        A[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        A0 = B0 / kb0
        A *= A0

    elif Flow_kind==3:
        # setting up flow (warning, this field has Vrms=sqrt(2))
        v[:, 0] = np.sin(kv0 * pos[:, 0])
        v[:, 1] = np.sin(kv0 * pos[:, 1])
        v[:, 2] = 2 * Vz_factor * np.cos(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
        v *= V0

        B[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        B[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        B[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        B *= B0

        A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        A[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        A[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        A0 = B0 / kb0
        A *= A0

    elif Flow_kind==4:
        # setting up flow (warning, this field has Vrms=sqrt(2))
        v[:, 0] = np.sin(kv0 * pos[:, 0])
        v[:, 1] = np.sin(kv0 * pos[:, 1])
        v[:, 2] = Vz_factor * np.sin(kv0 * (pos[:, 0] + pos[:, 1]))
        v *= V0

        B[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        B[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        B[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        B *= B0

        A[:, 0] = np.sin(kb0 * pos[:, 2]) + np.cos(kb0 * pos[:, 1])
        A[:, 1] = np.sin(kb0 * pos[:, 0]) + np.cos(kb0 * pos[:, 2])
        A[:, 2] = np.sin(kb0 * pos[:, 1]) + np.cos(kb0 * pos[:, 0])
        A0 = B0 / kb0
        A *= A0
    else:
        print('Wrong Flow kind. Use values 0-5')
        Continue = False
    ###---------------------------###
else:
    from swiftsimio import load
    from swiftsimio.visualisation.slice import slice_gas
    from swiftsimio.visualisation.rotation import rotation_matrix_from_vector

    vol = L ** 3
    #Pue path to IC snapshots here
    filename = '../ICfiles/g32_randB_withVP.hdf5'
    data = load(filename)
    #print(data.metadata.gas_properties.field_names)

    pos = data.gas.coordinates[:].value
    rho = data.gas.densities.value
    h = data.gas.smoothing_lengths.value
    v = data.gas.velocities.value
    m = data.gas.masses.value
    P = data.gas.pressures.value
    B = data.gas.magnetic_flux_densities.value
    A = data.gas.magnetic_vector_potentials.value
    u = data.gas.internal_energies.value
    ids = data.gas.particle_ids.value
    N=len(h)


if Continue:
    # File
    print('Wrinting ICs ...')
    fileOutput = h5py.File(fileOutputName, "w")

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
