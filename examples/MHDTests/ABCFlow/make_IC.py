#############################
# This file is part of SWIFT
# Copyright
#############################

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Parameters
rho0 = 1.0
cs2 = 3025.0
gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))
Bi_fraction = 1e-8

# output file
fileOutputName = "ABCFlow.hdf5"

############################################################### Random B and A vector field generator
def generate_random_vectors(N_of_vectors, randomize_length=True):
    RVF = []
    for i in range(N_of_vectors):
        phi = np.random.rand() * 2 * np.pi
        theta = np.random.rand() * np.pi
        r = np.array(
            [np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)]
        )
        if randomize_length:
            r *= np.random.rand()
        RVF += [r]
    return np.array(RVF)


def generate_k(Npart, Lbox, min_sine_res=6):
    # calculate maximal wavenumber
    nmax = int(Npart ** (1 / 3) / min_sine_res)
    kmax = 2 * np.pi / Lbox * nmax
    k = []
    for nx in range(-nmax, nmax + 1):
        for ny in range(-nmax, nmax + 1):
            for nz in range(-nmax, nmax + 1):
                if np.sqrt(nx ** 2 + ny ** 2 + nz ** 2) <= nmax:
                    k.append([nx, ny, nz])
    k = np.array(k)
    k = k.astype("float64")
    k *= 2 * np.pi / Lbox
    return k, kmax


def generate_random_complex_phases(N_of_phases, generate_vector_phase=True):
    if generate_vector_phase:
        Phx = np.random.rand(N_of_phases) * 2 * np.pi
        Phy = np.random.rand(N_of_phases) * 2 * np.pi
        Phz = np.random.rand(N_of_phases) * 2 * np.pi
        phase_array = np.array([[Phx[i], Phy[i], Phz[i]] for i in range(N_of_phases)])
        z = np.cos(phase_array) + 1.0j * np.sin(phase_array)
    else:
        Ph = np.random.rand(N_of_phases) * 2 * np.pi
        z = np.cos(Ph) + 1.0j * np.sin(Ph)
    return z


def generate_Ak(kvec, kmax):
    a = generate_random_vectors(len(kvec), randomize_length=False)
    b = np.cross(a, kvec)
    kabs = np.linalg.norm(kvec, axis=1)
    babs = np.linalg.norm(b, axis=1)

    for i in range(len(b)):
        if babs[i] != 0:
            b[i] /= babs[i]
        else:
            b[i] = a[i]

    Cph = generate_random_complex_phases(len(kvec))
    amps = np.random.rand(len(kvec))
    Ak = np.zeros(kvec.shape) + 1.0j * np.zeros(kvec.shape)

    for i in range(len(kvec)):
        Ak[i] = Cph[i] * b[i] * amps[i]

    return Ak


def find_Magnetic_Field(pos, kvec, Ak):
    A = np.zeros(pos.shape)
    B = np.zeros(pos.shape)
    Ph = 2 * np.pi * np.random.rand(len(kvec))
    for i in range(len(pos)):
        for j in range(3):
            A[i, j] = np.real(
                np.sum(
                    Ak[:, j]
                    * (
                        np.exp(
                            1.0j
                            * (
                                kvec[:, 0] * pos[i, 0]
                                + kvec[:, 1] * pos[i, 1]
                                + kvec[:, 2] * pos[i, 2]
                                + Ph
                            )
                        )
                    )[:, None]
                )
            )
            B[i, j] = np.real(
                1.0j
                * np.sum(
                    np.cross(kvec[:], Ak[:])[:, j]
                    * (
                        np.exp(
                            1.0j
                            * (
                                kvec[:, 0] * pos[i, 0]
                                + kvec[:, 1] * pos[i, 1]
                                + kvec[:, 2] * pos[i, 2]
                                + Ph
                            )
                        )
                    )[:, None]
                )
            )
        print(
            "Generating magnetic field: "
            + str(int((i + 1) / len(pos) * 100))
            + "% complete",
            end="\r",
        )
    return A, B


def normalize_Magnetic_Field(A, B, B0):
    norms = np.linalg.norm(B, axis=1)
    normalization_factor = B0 / np.sqrt(np.sum(norms ** 2))
    return A * normalization_factor, B * normalization_factor


############################################################################


# smoothing kernel optimized for numpy Price 1012.1885 (6)
def open_IAfile(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLength"][:]
    return pos, h


def add_other_particle_properties(
    pos, h, a, b, c, V0, kb, kv, Vz_factor, L, field_type
):
    if field_type != "load_from_file":
        vol = L ** 3
        N = len(h)

        # initializing arrays with particle properties
        v = np.zeros((N, 3))
        B = np.zeros((N, 3))
        A = np.zeros((N, 3))
        ids = np.linspace(1, N, N)
        m = np.ones(N) * rho0 * vol / N
        u = np.ones(N) * u0

        # setting constants
        kv0 = 2 * np.pi / L * kv
        kb0 = 2 * np.pi / L * kb
        Beq0 = np.sqrt(rho0) * V0
        B0 = Bi_fraction * Beq0
        Norm = 1 / np.sqrt(a ** 2 + b ** 2 + c ** 2)

        # rescaling the box to size L
        pos *= L

        # setting the velocity profile
        v[:, 0] = a * np.sin(kv0 * pos[:, 2]) + c * np.cos(kv0 * pos[:, 1])
        v[:, 1] = b * np.sin(kv0 * pos[:, 0]) + a * np.cos(kv0 * pos[:, 2])
        v[:, 2] = c * np.sin(kv0 * pos[:, 1]) + b * np.cos(kv0 * pos[:, 0])
        v[:, 2] *= Vz_factor
        v *= V0 * Norm

        # setting the initial magnetic field

        if field_type == "one_mode":
            B[:, 0] = -(np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1]))
            B[:, 1] = -(np.cos(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2]))
            B[:, 2] = -(np.cos(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0]))
            B *= B0

            A[:, 0] = np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1])
            A[:, 1] = np.sin(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2])
            A[:, 2] = np.sin(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0])
            A0 = B0 / kb0
            A *= A0
        elif field_type == "random":
            # kvec,kmax = generate_k(N,L)
            # Ak = generate_Ak(kvec,kmax)
            # A,B = find_Magnetic_Field(pos,kvec,Ak)
            B = generate_random_vectors(N, randomize_length=True)
            A = B * L / (2 * np.pi)
            A, B = normalize_Magnetic_Field(A, B, B0)
        else:
            print("Error: wrong field type. Should be one_mode or random")
    else:

        vol = L ** 3
        # Put path to IC snapshots here
        filename = "./ICfiles/g32_randB_withVP.hdf5"
        # read the variables of interest from the snapshot file
        pos = None
        h = None
        v = None
        B = None
        A = None
        ids = None
        with h5py.File(filename, "r") as handle:
            print(handle["PartType0"].keys())
            pos = handle["/PartType0/Coordinates"][:]
            Lbox = handle["Header"].attrs.get("BoxSize")
            h = handle["/PartType0/SmoothingLengths"][:]
            v = handle["PartType0/Velocities"][:]
            m = handle["PartType0/Masses"][:]
            u = handle["PartType0/InternalEnergies"][:]
            B = handle["PartType0/MagneticFluxDensities"][:]
            A = handle["PartType0/MagneticVectorPotentials"][:]
            ids = handle["PartType0/ParticleIDs"][:]

        pos *= L / Lbox[0]
        h *= L / Lbox[0]
        V0_from_snap = np.sqrt(np.mean(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2))
        v *= V0 / V0_from_snap
        Beq0 = np.sqrt(rho0) * V0
        B0 = Bi_fraction * Beq0
        A, B = normalize_Magnetic_Field(A, B, B0)
        # data = load(filename)
        # print(data.metadata.gas_properties.field_names)

        # pos = data.gas.coordinates[:].value
        # rho = data.gas.densities.value
        # h = data.gas.smoothing_lengths.value
        # v = data.gas.velocities.value
        # m = data.gas.masses.value
        # P = data.gas.pressures.value
        # B = data.gas.magnetic_flux_densities.value
        # A = data.gas.magnetic_vector_potentials.value
        # u = data.gas.internal_energies.value
        # ids = data.gas.particle_ids.value
        # N = len(h)

    return pos, h, v, B, A, ids, m, u


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
        default="./IAfiles/glassCube_32.hdf5",
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
        default=2 * np.pi,
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
        "-A", "--A", help="constant A of the ABC flow", default=1.0, type=float
    )
    parser.add_argument(
        "-B", "--B", help="constant B of the ABC flow", default=1.0, type=float
    )
    parser.add_argument(
        "-C", "--C", help="constant C of the ABC flow", default=1.0, type=float
    )
    parser.add_argument(
        "-ft",
        "--field_type",
        help="How to generate a field: one_mode, several_modes or random",
        default="random",  #'load_from_file',#'random',
        type=str,
    )
    args = parser.parse_args()
    pos, h = open_IAfile(args.IA_path)
    pos, h, v, B, A, ids, m, u = add_other_particle_properties(
        pos,
        h,
        args.A,
        args.B,
        args.C,
        args.rms_velocity,
        args.magnetic_wavevector,
        args.velocity_wavevector,
        args.Vz_factor,
        args.boxsize,
        args.field_type,
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
