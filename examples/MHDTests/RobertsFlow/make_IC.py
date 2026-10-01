#############################
# This file is part of SWIFT
# Copyright


# Note: the code below contains AI generated functions:
# spectrum_equal_mode_amplitude
# spectrum_powerlaw
# generate_random_magnetic_field_from_spectra

#############################

import h5py
import numpy as np
#import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

# Parameters
rho0 = 1.0
cs2 = 1.0 #3025.0
gamma = 5.0 / 3.0
u0 = cs2 / (gamma * (gamma - 1))
Bi_fraction = 1e-6

# output file
fileOutputName = "RobertsFlow.hdf5"

def spectrum_powerlaw(k, n, kmin, kmax):
    """
    Spectral amplitude shape.
    The overall normalization does not matter.
    """

    if n!=0:
        P = np.zeros_like(k)
        mask = ( ( k >= kmin) & (k <= kmax) )
        P[mask] = k[mask]**n
    else:
        P = np.zeros_like(k)
        mask = ( (k >= kmin) & (k <= kmax))
        P[mask] = 1.0

    return P


def generate_random_magnetic_field_from_spectra(
    positions,
    boxsize,
    Ngrid,
    Brms,
    spectrum,
    seed=None,
    kmin=None,
    kmax=None,
    return_grid=False,
):
    """
    Generate a divergence-free random magnetic field
    with a prescribed spectrum.

    Parameters
    ----------
    positions : ndarray, shape (Nparticles, 3)
        Particle coordinates.

    boxsize : array-like, shape (3,)
        Periodic box dimensions [Lx, Ly, Lz].

    Ngrid : int or tuple of 3 ints
        Number of grid cells per dimension.

    Brms : float
        Desired RMS magnetic field.

    spectrum : callable
        Function P(k), specifying the spectral shape.

        The normalization of P(k) is arbitrary because
        the final field is normalized to Brms.

    seed : int, optional
        Random seed.

    kmin, kmax : float, optional
        Minimum and maximum wavenumbers.
        Modes outside this interval are removed.

    return_grid : bool
        If True, also return the grid field.

    Returns
    -------
    B_particles : ndarray, shape (Nparticles, 3)
        Magnetic field at particle positions.

    B_grid : ndarray, shape (Nx, Ny, Nz, 3), optional
        Magnetic field on the grid.

    """

    rng = np.random.default_rng(seed)

    positions = np.asarray(positions, dtype=float)
    boxsize = np.asarray(boxsize, dtype=float)

    if np.isscalar(Ngrid):
        Ngrid = (Ngrid, Ngrid, Ngrid)

    Nx, Ny, Nz = Ngrid
    Lx, Ly, Lz = boxsize

    dx = boxsize / np.array(Ngrid)

    # --------------------------------------------------
    # 1. Fourier wavevectors
    # --------------------------------------------------

    kx = 2.0 * np.pi * np.fft.fftfreq(Nx, d=dx[0])
    ky = 2.0 * np.pi * np.fft.fftfreq(Ny, d=dx[1])
    kz = 2.0 * np.pi * np.fft.fftfreq(Nz, d=dx[2])

    KX, KY, KZ = np.meshgrid(
        kx, ky, kz, indexing="ij"
    )

    K2 = KX**2 + KY**2 + KZ**2
    K = np.sqrt(K2)

    # --------------------------------------------------
    # 2. Generate random vector field in Fourier space
    # --------------------------------------------------

    # Independent real Gaussian random fields.
    # The FFT of a real field has Hermitian symmetry.
    noise = rng.normal(
        size=(Nx, Ny, Nz, 3)
    )

    Bk = np.empty(
        (Nx, Ny, Nz, 3),
        dtype=complex
    )

    for component in range(3):
        Bk[..., component] = np.fft.fftn(
            noise[..., component]
        )

    # --------------------------------------------------
    # 3. Apply desired spectral amplitude
    # --------------------------------------------------

    Pk = np.asarray(spectrum(K), dtype=float)

    Pk[~np.isfinite(Pk)] = 0.0
    Pk[Pk < 0] = 0.0

    if kmin is not None:
        Pk[K < kmin] = 0.0

    if kmax is not None:
        Pk[K > kmax] = 0.0

    # No magnetic field at k=0.
    Pk[K == 0] = 0.0

    Bk *= np.sqrt(Pk)[..., None]

    # --------------------------------------------------
    # 4. Project Fourier modes perpendicular to k
    # --------------------------------------------------

    # Projection:
    #
    # B_perp = B - k (k dot B) / |k|^2
    #
    # This guarantees k dot B_perp = 0.

    dot = (
        KX * Bk[..., 0]
        + KY * Bk[..., 1]
        + KZ * Bk[..., 2]
    )

    Bk[..., 0] -= KX * dot / np.where(K2 == 0, 1.0, K2)
    Bk[..., 1] -= KY * dot / np.where(K2 == 0, 1.0, K2)
    Bk[..., 2] -= KZ * dot / np.where(K2 == 0, 1.0, K2)

    Bk[K2 == 0] = 0.0


    # --------------------------------------------------
    # 5. Construct vector potential in Coulomb gauge
    # --------------------------------------------------
    
    # B_k = i k x A_k
    #
    # Therefore:
    #
    # A_k = i (k x B_k) / k^2
    #
    # This gives:
    #
    # k . A_k = 0
    #
    # i k x A_k = B_k
    
    K_cross_B_x = KY * Bk[..., 2] - KZ * Bk[..., 1]
    K_cross_B_y = KZ * Bk[..., 0] - KX * Bk[..., 2]
    K_cross_B_z = KX * Bk[..., 1] - KY * Bk[..., 0]
    
    Ak = np.empty_like(Bk)
    
    K2_safe = np.where(K2 == 0.0, 1.0, K2)
    
    Ak[..., 0] = 1j * K_cross_B_x / K2_safe
    Ak[..., 1] = 1j * K_cross_B_y / K2_safe
    Ak[..., 2] = 1j * K_cross_B_z / K2_safe
    
    # No k=0 vector-potential mode
    Ak[K2 == 0] = 0.0

    # --------------------------------------------------
    # 6. Transform to real space
    # --------------------------------------------------

    B_grid = np.empty(
        (Nx, Ny, Nz, 3),
        dtype=float
    )
    A_grid = np.empty(
        (Nx, Ny, Nz, 3),
        dtype=float
    )

    for component in range(3):
        B_grid[..., component] = np.fft.ifftn(
            Bk[..., component]
        ).real
        A_grid[..., component] = np.fft.ifftn(
            Ak[..., component]
        ).real


    # --------------------------------------------------
    # 7. Normalize to desired Brms
    # --------------------------------------------------

    current_Brms = np.sqrt(
        np.mean(np.sum(B_grid**2, axis=-1))
    )

    if current_Brms == 0.0:
        raise ValueError(
            "Generated field has zero RMS. "
            "Check the spectrum and k limits."
        )

    B_grid *= Brms / current_Brms
    A_grid *= Brms / current_Brms

    # --------------------------------------------------
    # 8. Interpolate B field to SPH particles
    # --------------------------------------------------

    grid_axes = [
        np.arange(Nx) * dx[0],
        np.arange(Ny) * dx[1],
        np.arange(Nz) * dx[2],
    ]

    # Periodic wrapping
    positions_wrapped = np.mod(
        positions, boxsize
    )

    B_particles = np.empty_like(positions)

    for component in range(3):

        interpolator = RegularGridInterpolator(
            grid_axes,
            B_grid[..., component],
            method="linear",
            bounds_error=False,
            fill_value=None,
        )

        # Need periodic interpolation at the boundary.
        # Extend one grid point by periodic wrapping.
        field = B_grid[..., component]

        field_ext = np.pad(
            field,
            ((0, 1), (0, 1), (0, 1)),
            mode="wrap"
        )

        axes_ext = [
            np.append(grid_axes[0], boxsize[0]),
            np.append(grid_axes[1], boxsize[1]),
            np.append(grid_axes[2], boxsize[2]),
        ]

        interpolator = RegularGridInterpolator(
            axes_ext,
            field_ext,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )

        B_particles[:, component] = interpolator(
            positions_wrapped
        )

    # --------------------------------------------------
    # 9. Interpolate A_grid to SPH particles
    # --------------------------------------------------
    
    grid_axes = [
        np.arange(Nx) * dx[0],
        np.arange(Ny) * dx[1],
        np.arange(Nz) * dx[2],
    ]
    
    positions_wrapped = np.mod(
        positions,
        boxsize
    )
    
    A_particles = np.empty_like(positions)
    
    for component in range(3):
    
        field = A_grid[..., component]
    
        # Periodic extension
        field_ext = np.pad(
            field,
            ((0, 1), (0, 1), (0, 1)),
            mode="wrap"
        )
    
        axes_ext = [
            np.append(grid_axes[0], boxsize[0]),
            np.append(grid_axes[1], boxsize[1]),
            np.append(grid_axes[2], boxsize[2]),
        ]
    
        interpolator = RegularGridInterpolator(
            axes_ext,
            field_ext,
            method="linear",
            bounds_error=False,
            fill_value=None,
        )
    
        A_particles[:, component] = interpolator(
            positions_wrapped
        )

    print('Random magnetic field spectrum added')

    if return_grid:
        return B_particles, B_grid, A_particles, A_grid

    return B_particles, A_particles



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


def normalize_Magnetic_Field(A, B, B0):
    norms = np.linalg.norm(B, axis=1)
    normalization_factor = B0 / np.sqrt(np.sum(norms ** 2))
    return A * normalization_factor, B * normalization_factor


############################################################################


def open_IAfile(path_to_file):
    IAfile = h5py.File(path_to_file, "r")
    pos = IAfile["/PartType0/Coordinates"][:, :]
    h = IAfile["/PartType0/SmoothingLength"][:]
    return pos, h


def add_other_particle_properties(
    pos, h, V0, kb, kv, Vz_factor, L, field_type, Flow_kind
):
    # One can start from any .hdf5 outputs of the simulation 
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

        # rescaling the box to size L
        pos *= L

        # setting the velocity profile
        if Flow_kind == 1:
            v[:, 0] = np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
            v[:, 1] = -np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
            v[:, 2] = (
                Vz_factor
                * np.sqrt(2)
                * np.sin(kv0 * pos[:, 1])
                * np.sin(kv0 * pos[:, 0])
            )
        elif Flow_kind == 2:
            v[:, 0] = np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
            v[:, 1] = -np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
            v[:, 2] = (
                Vz_factor
                * np.sqrt(2)
                * np.cos(kv0 * pos[:, 0])
                * np.cos(kv0 * pos[:, 1])
            )
        elif Flow_kind == 3:
            v[:, 0] = np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
            v[:, 1] = -np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
            v[:, 2] = (
                Vz_factor
                * (np.cos(2 * kv0 * pos[:, 0]) + np.cos(2 * kv0 * pos[:, 1]))
                / np.sqrt(2)
            )
        elif Flow_kind == 4:
            v[:, 0] = np.sin(kv0 * pos[:, 0]) * np.cos(kv0 * pos[:, 1])
            v[:, 1] = -np.sin(kv0 * pos[:, 1]) * np.cos(kv0 * pos[:, 0])
            v[:, 2] = Vz_factor * np.sin(kv0 * pos[:, 0])
        else:
            print("Wrong Flow kind. Use values 1-4")

        v *= V0

        # setting the initial magnetic field

        if field_type == "one_mode":
            B[:, 0] = -(np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1]))
            B[:, 1] = -(np.sin(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2]))
            B[:, 2] = -(np.sin(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0]))
            B *= B0
           
            meanB = np.mean(B,axis=0)
            B -= meanB[None,:]

            A[:, 0] = np.sin(kb0 * pos[:, 2]) - np.cos(kb0 * pos[:, 1])
            A[:, 1] = np.sin(kb0 * pos[:, 0]) - np.cos(kb0 * pos[:, 2])
            A[:, 2] = np.sin(kb0 * pos[:, 1]) - np.cos(kb0 * pos[:, 0])
            A0 = B0 / kb0
            A *= A0
        if field_type == "BxSin":
            B[:, 0] = np.sin(kb0 * pos[:, 2])

            meanB = np.mean(B,axis=0)
            B -= meanB[None,:]

            B *= B0
            A[:, 1] = - np.cos(kb0 * pos[:, 2])
            A0 = B0 / kb0
            A *= A0
 
        elif field_type == "random":

            Npx = int((len(B))**(1/3))

            # generate magnetic field with k^-3 spectrum and a cutoff at the Nyquist scale, with d_res = kernel cutoff radius ~ 3.2 interparticle spacing
            nB = -3.0

            B, B_grid = generate_random_magnetic_field_from_spectra(
                positions=pos,
                boxsize=np.array([L,L,L]),
                Ngrid=4*Npx,
                Brms=B0,
                spectrum=lambda k: spectrum_powerlaw(
                    k,
                    n=-nB,
                    kmin = 2.0 * np.pi / L,
                    kmax = np.pi / L * (Npx / 3.2),
                ),
                seed=1234,
                return_grid=True,
                )

            A = B * L / (2 * np.pi)
        else:
            print("Error: wrong field type. Should be one_mode or random")

    # One can start from any .hdf5 outputs of the simulation 
    else:
        vol = L ** 3
        # Put path to IC snapshots here
        filename ="./ICfiles/rf2d_o1_npar2_g64_randB_withVP.hdf5" #"./ICfiles/rf4d_o1_g64_randB_withVP.hdf5" #"./ICfiles/RF2_pattern_0090.hdf5" #"./ICfiles/rf2d_o1_g64_randB_withVP.hdf5" #"./ICfiles/rf2d_g128_randB_withVP.hdf5"
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
            try:
                A = handle["PartType0/MagneticVectorPotentials"][:]
            except:
                print('No VP, continue')
                A = np.zeros(B.shape)
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

def stack_boxes(pos,h,npar,nper):
    h_stack = []
    pos_stack = []
    for i in range(npar):
        for j in range(npar):
            for k in range(nper):
                h_stack.append(h[:])
                pos_stack.append(pos[:]+np.array([i,j,k])[None,:]) 

    pos_stack = np.concatenate(pos_stack)
    h_stack = np.concatenate(h_stack)

   # N = len(h)

    #cx, cy, cz = npar, npar, nper

    #pos_stack = np.zeros((int(N * cx * cy * cz), 3))
   # h_stack = np.zeros(int(N * cx * cy * cz))
   # N_stack = N * cx * cy * cz

    #c0 = 0
    #c1 = N
    #for i in range(cx):
    #    for j in range(cy):
    #        for k in range(cz):
    #            pos_stack[c0:c1, 0] = pos[:, 0] + i
    #            pos_stack[c0:c1, 1] = pos[:, 1] + j
    #            pos_stack[c0:c1, 2] = pos[:, 2] + k
    #            h_stack[c0:c1] = h[:]
    #            c0 += N
    #            c1 += N

    return pos_stack, h_stack


def deform_boxes(pos,h,Lparmul,Lpermul):
    pos[:,0]*=Lparmul
    pos[:,1]*=Lparmul
    pos[:,2]*=Lpermul
    h[:]*= np.sqrt(2*Lparmul**2+Lpermul**2)
    return pos,h 

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
        default="./IAfiles/glassCube_16.hdf5",
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
        "-ft",
        "--field_type",
        help="How to generate a field: one_mode, several_modes or random",
        default="random", #"BxSin", #"one_mode",#"load_from_file", #"BxSin", #"one_mode", #"load_from_file",  #'load_from_file',#'random',
        type=str,
    )
    parser.add_argument(
        "-fk",
        "--flow_kind",
        help="RobertsFlow has four flow kinds: 1-4",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-npar",
        "--npar_box",
        help="number of boxes in xy plane",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-nper",
        "--nper_box",
        help="number of boxes in z direction",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-lparmul",
        "--lparmultiplier",
        help="multiply boxsize in xy plane by this amount",
        default=1,
        type=float,
    )
    parser.add_argument(
        "-lpermul",
        "--lpermultiplier",
        help="number of boxes in z plane by this amount",
        default=1,
        type=float,
    )

    parser.add_argument(
        "-vtk",
        "--to_vtk",
        help="wether to save result to .vtk file, 1 or 0",
        default=0,
        type=int,
    )



    args = parser.parse_args()
    pos, h = open_IAfile(args.IA_path)

    pos, h = stack_boxes(
        pos,
        h,
        args.npar_box,
        args.nper_box,
    )

    print(pos, h)

    pos, h = deform_boxes(
        pos,
        h,
        args.lparmultiplier,
        args.lpermultiplier
    )

    # Find unique rows and their counts
    unique_pos, counts = np.unique(pos, axis=0, return_counts=True)

    # Filter positions that appear more than once
    repeating_pos = unique_pos[counts > 1]

    print("Repeating positions:")
    print(repeating_pos)

    pos, h, v, B, A, ids, m, u = add_other_particle_properties(
        pos,
        h,
        args.rms_velocity,
        args.magnetic_wavevector,
        args.velocity_wavevector,
        args.Vz_factor,
        args.boxsize,
        args.field_type,
        args.flow_kind,
    )

    L = args.boxsize
    cpar = args.npar_box*args.lparmultiplier
    cper = args.nper_box*args.lpermultiplier
    Lbox = [L*cpar, L*cpar, L*cper]
    N = len(h)

    # try to separate the modes
    #Brms = np.sqrt(np.mean(np.linalg.norm(B,axis=1)**2))

    #B[:,0] = 0.0
    #B[:,1] = 0.0
    #B[:,2] = 0.0

    # excite Bx
    #B[:,0] += Brms * np.sin(2 * np.pi * pos[:,2]/Lbox[2])

    # excite By
    #B[:,1] += Brms * np.sin(2 * np.pi * pos[:,2]/Lbox[2])

    # File
    try:
        fileOutput = h5py.File(fileOutputName, "w")
        print("Wrinting ICs ...")
        # Header
        grp = fileOutput.create_group("/Header")
        grp.attrs["BoxSize"] = Lbox  #####
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

    # Save particle arrangement to .vtk file to open with paraview
    if args.to_vtk==1:
        import vtk

        vtk_points = vtk.vtkPoints()
        for x, y, z in pos:
            vtk_points.InsertNextPoint(x, y, z)
        poly_data = vtk.vtkPolyData()
        poly_data.SetPoints(vtk_points)
        writer = vtk.vtkPolyDataWriter()
        writer.SetFileName("output.vtk")
        writer.SetInputData(poly_data)
        writer.Write()  # one can use paraview to watch how particles are arranged


