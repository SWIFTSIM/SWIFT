import os
import re
import numpy as np

# import h5py
import matplotlib.pyplot as plt
import argparse
import scipy as sp

from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

# SWIFT-only dependencies: Pencil snapshots can be processed without them
try:
    import unyt
    from swiftsimio import load
    from swiftsimio.visualisation.slice import slice_gas
except ImportError as _e:
    _swift_import_error = _e
else:
    _swift_import_error = None

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",  # Computer Modern
    }
)


# Parse command line arguments
# argparser = argparse.ArgumentParser()
# argparser.add_argument("input")
# argparser.add_argument("output")
# argparser.add_argument("resolution")
# args = argparser.parse_args()


# Load snapshot
# filename = args.input


def periodic_continue(Q):
    """
    Periodically continue the image
    Parameters: Q: quantity, with size (N_x, N_y)
    Returns:
                Qp: quantity, with size (N_x+1, N_y+1)
    """

    M, N = Q.shape
    # Q = np.arange(M*N).reshape(M, N)  # example data

    # Periodic continuation
    Qp = np.zeros((M + 1, N + 1), dtype=Q.dtype)
    Qp[:-1, :-1] = Q  # copy original
    Qp[-1, :-1] = Q[0, :]  # wrap last row
    Qp[:-1, -1] = Q[:, 0]  # wrap last column
    Qp[-1, -1] = Q[0, 0]  # wrap bottom-right corner
    return Qp


def make_2D_slice_z(data, height, res_xy):
    """
    Making SPH B slice of the box at specific height
    Parameters: data: SPH snapshot data
                height: slice height
    Returns:
                dict: B - magnetic field array of size (N_x, N_y, 3), z - slice height
    """

    Lbox = data.metadata.boxsize

    # compensate half-grid shift
    dx = Lbox[0] / (2 * res_xy)
    dy = Lbox[1] / (2 * res_xy)

    # select region
    visualise_region_xy = [
        0 * unyt.cm - dx,
        Lbox[0] - dx,
        0 * unyt.cm - dy,
        Lbox[1] - dy,
    ]

    common_arguments_xy = dict(
        data=data,
        resolution=res_xy,
        parallel=True,
        region=visualise_region_xy,
        z_slice=height,
        periodic=True,
    )

    mass_map_ij = slice_gas(**common_arguments_xy, project="masses")
    mass_weighted_Bx_map_ij = slice_gas(
        **common_arguments_xy, project="mass_weighted_Bx"
    )
    mass_weighted_By_map_ij = slice_gas(
        **common_arguments_xy, project="mass_weighted_By"
    )
    mass_weighted_Bz_map_ij = slice_gas(
        **common_arguments_xy, project="mass_weighted_Bz"
    )
    Bx_map_ij = mass_weighted_Bx_map_ij / mass_map_ij
    By_map_ij = mass_weighted_By_map_ij / mass_map_ij
    Bz_map_ij = mass_weighted_Bz_map_ij / mass_map_ij

    # Invert y axis to convert ij map to xy map
    Bx_map_xy = Bx_map_ij  # np.flipud(Bx_map_ij)
    By_map_xy = By_map_ij  # np.flipud(By_map_ij)
    Bz_map_xy = Bz_map_ij  # np.flipud(Bz_map_ij)

    B = np.stack((Bx_map_xy, By_map_xy, Bz_map_xy), axis=-1)
    return {"B": B, "z": height}


def make_SPH_grid_B(data, res_xy, res_z):
    """
    Rendering B field from SPH particle disrtibution
    Parameters: data: SPH snapshot data
                res_xy: number of samples along x or y
                res_z: number of samples along x or y
    Returns:
                Bx,By,Bz: 3 3d arrays of shape (N_x, N_y, N_z)
    """

    # Retrieve particle attributes of interest
    B = data.gas.magnetic_flux_densities
    Lbox = data.metadata.boxsize
    # Generate mass weighted maps of quantities of interest
    data.gas.mass_weighted_Bx = data.gas.masses * B[:, 0]
    data.gas.mass_weighted_By = data.gas.masses * B[:, 1]
    data.gas.mass_weighted_Bz = data.gas.masses * B[:, 2]

    slices = np.linspace(0, Lbox[2].value, res_z, endpoint=False) * Lbox[2].units

    Bx = []
    By = []
    Bz = []

    for slice in slices:
        Bx.append(make_2D_slice_z(data, slice, res_xy)["B"][:, :, 0])
        By.append(make_2D_slice_z(data, slice, res_xy)["B"][:, :, 1])
        Bz.append(make_2D_slice_z(data, slice, res_xy)["B"][:, :, 2])

    # make with shape N_x, N_y, N_z
    Bx = np.array(Bx).transpose(1, 2, 0)
    By = np.array(By).transpose(1, 2, 0)
    Bz = np.array(Bz).transpose(1, 2, 0)

    return Bx, By, Bz


def compute_fourier_transform_1d(Q, dx):
    """
    Applying fast Fourier transformation to a quantity Q
    Parameters: Q: 1D array
                dx: spacial discretization
    Returns:
                k: wavenumbers for each mode
                Q_k: Fourier image of Q
    """

    # Grid size
    N = len(Q)

    # Compute Fourier transforms
    # Q_k = np.fft.fftn(Q) * dx / Lbox[2]
    Q_k = np.fft.fft(Q) * dx / Lbox[2]

    # Compute the corresponding wavenumbers
    k = np.fft.fftfreq(N, d=dx) * 2 * np.pi

    return k, Q_k


def transform_B_along_z(
    Bx,
    By,
    Bz,
    res_xy,
    res_z,
):
    """
    Applying Fourier transformation along z to the B field array
    Parameters: Bx,By,Bz: 3D arrays with B field components
                res_xy: number of samples along x or y
                res_z: number of samples along z
    Returns:
                k: wavenumbers for each mode
                b_xy: array of shape (N_x,N_y, N_modes, N_dim=3) - Fourier image of B
    """

    dz = Lbox[2] / res_z

    b_xy = []
    for i in range(res_xy):
        b_xy_row = []
        for j in range(res_xy):
            k, bx_k = compute_fourier_transform_1d(Bx[i, j, :], dz)
            k, by_k = compute_fourier_transform_1d(By[i, j, :], dz)
            k, bz_k = compute_fourier_transform_1d(Bz[i, j, :], dz)
            b_k = np.stack((bx_k, by_k, bz_k), axis=-1)
            b_xy_row.append(b_k)
        b_xy.append(b_xy_row)

    return k, np.array(b_xy)


def save_to_sav(savedat: dict, outputfilename: str = "b0perp.sav"):
    """
    Save the data from dict to .sav file
    Parameters: savedat: data to save, dict of the form {'item1':array, ...}
                outputfilename: output filename
    Returns:
    """

    # Save in IDL .sav format with the variable name 'b0perp'
    sp.io.savemat(outputfilename, savedat)


def load_cmap(img_addr: str, cmapname: str):
    """
    Load custom colormaps (loads an image and creates colormap).
    Parameters: img_addr: adress of the colormap image
                cmapname: name of the colormap
    Returns: colormap object
    """

    img = imread(img_addr)

    colors_from_img = img[0, :, :]

    size = len(colors_from_img)

    res = LinearSegmentedColormap.from_list(cmapname, colors_from_img, N=size)
    return res


def generate_test_field(data):
    """
    Generate a test field for testing purposes.
    Parameters: data: snapshot data
    Returns: test_field: a 3D array
    """

    import numpy as np

    # res = 64 # 16 #64 #16
    Nz = res

    # Define box size (to mimic Lbox from snapshot metadata)
    Lbox = data.metadata.boxsize

    # np.array([1.0, 1.0, 1.0])  # in arbitrary units

    # Create coordinates
    x = np.linspace(0, Lbox[0].value, res, endpoint=False) * Lbox.units
    y = np.linspace(0, Lbox[1].value, res, endpoint=False) * Lbox.units
    z = np.linspace(0, Lbox[2].value, Nz, endpoint=False) * Lbox.units

    Z, X, Y = np.meshgrid(z, x, y, indexing="ij")

    # Example magnetic fields
    Bx = np.zeros_like(Y)  # uniform field along x
    By = np.zeros_like(
        Y
    )  # np.zeros_like(Y)       #1.0*np.cos(2*np.pi*Z/Lbox[2])               # no y-component
    Bz = np.cos(
        2 * np.pi * Z / Lbox[2] + X.value
    )  # np.zeros_like(Z)                              #np.cos(2*np.pi*Z/Lbox[2])*np.cos(2*np.pi*X/Lbox[2]) #np.zeros_like(Y) #np.cos(2*np.pi*Z/Lbox[2]+0.5*np.pi*np.cos(2*np.pi*X/Lbox[0]))  #np.zeros_like(Z) #1.0*np.cos(2*np.pi*Z/Lbox[2])               # np.zeros_like(Z)
    B_data = (Bx, By, Bz) * data.gas.magnetic_flux_densities.units
    print(Bx)
    return B_data


# def read_fits_B(input_filename,
#    Lbox = [2*np.pi,2*np.pi, 4*np.pi],
#    ):
#    hdul = fits.open(input_filename)
#    print(hdul.info())
#    primary_hdu = hdul[0]
#    data = primary_hdu.data
#    header = primary_hdu.header
#
#    Bx,By,Bz = data
#
#    # make with shape N_x, N_y, N_z
#    Bx = np.array(Bx).transpose(2, 1, 0)
#    By = np.array(By).transpose(2, 1, 0)
#    Bz = np.array(Bz).transpose(2, 1, 0)
#
#    return Bx,By,Bz, Lbox


def read_hdf5_B(input_filename, res_xy, res_z):
    if _swift_import_error is not None:
        raise ImportError(f"SWIFT snapshots need swiftsimio and unyt: {_swift_import_error}")
    # Load snapshot
    data = load(input_filename)

    time = np.round(data.metadata.time, 2)
    print("Plotting at ", time)

    Lbox = data.metadata.boxsize

    uu = data.gas.velocities
    vrms = np.sqrt(np.mean(uu[:, 0] ** 2 + uu[:, 1] ** 2 + uu[:, 2] ** 2))
    print("vrms", vrms)

    Bx, By, Bz = make_SPH_grid_B(data, res_xy, res_z)

    return Bx, By, Bz, Lbox.value


def _read_dim(filename):
    """Pencil dim.dat -> (mx, my, mz, precision, nghost (3,), last line (3,)).
    The last line is (nprocx, nprocy, nprocz) in data/dim.dat and (ipx, ipy, ipz) in procN/dim.dat."""
    with open(filename) as f:
        lines = [l.split() for l in f if l.strip()]
    mx, my, mz = map(int, lines[0][:3])
    return mx, my, mz, lines[1][0], tuple(map(int, lines[2][:3])), tuple(map(int, lines[3][:3]))


def read_pencil_var(datadir, varfile):
    """Vector potential from a Pencil snapshot data/proc*/varfile, read without the pencil package
    (fallback for read_VAR_B).

    Parameters: datadir: Pencil data directory (dim.dat, index.pro, proc0, ...)
                varfile: snapshot name, e.g. 'VAR2' or 'var.dat'
    Returns:
                A: array (3, nx, ny, nz) indexed [component, x, y, z], ghost zones removed
                t: snapshot time;  x, y, z: grid coordinates
    """
    mx, my, mz, prec, (gx, gy, gz), (npx, npy, npz) = _read_dim(os.path.join(datadir, "dim.dat"))
    nx, ny, nz = mx - 2 * gx, my - 2 * gy, mz - 2 * gz
    dtype = np.dtype("<f8" if prec == "D" else "<f4")
    with open(os.path.join(datadir, "index.pro")) as f:
        iax = int(re.search(r"^\s*iax\s*=\s*(\d+)", f.read(), re.M).group(1)) - 1
    A = np.zeros((3, nx, ny, nz))
    x, y, z = np.zeros(nx), np.zeros(ny), np.zeros(nz)
    for p in range(npx * npy * npz):
        pdir = os.path.join(datadir, f"proc{p}")
        lmx, lmy, lmz, _, _, (ipx, ipy, ipz) = _read_dim(os.path.join(pdir, "dim.dat"))
        lnx, lny, lnz = lmx - 2 * gx, lmy - 2 * gy, lmz - 2 * gz
        with open(os.path.join(pdir, varfile), "rb") as f:
            raw = f.read()
        # record 1: f(mx, my, mz, nvar), Fortran order
        n1 = int(np.frombuffer(raw[:4], "<i4")[0])
        nvar = n1 // (lmx * lmy * lmz * dtype.itemsize)
        fa = np.frombuffer(raw[4 : 4 + n1], dtype).reshape(nvar, lmz, lmy, lmx)
        # record 2: t, x(mx), y(my), z(mz), dx, dy, dz
        o = 4 + n1 + 4
        n2 = int(np.frombuffer(raw[o : o + 4], "<i4")[0])
        r2 = np.frombuffer(raw[o + 4 : o + 4 + n2], dtype)
        t = float(r2[0])
        i0, j0, k0 = ipx * lnx, ipy * lny, ipz * lnz
        a = fa[iax : iax + 3, gz : gz + lnz, gy : gy + lny, gx : gx + lnx]
        A[:, i0 : i0 + lnx, j0 : j0 + lny, k0 : k0 + lnz] = a.transpose(0, 3, 2, 1)
        x[i0 : i0 + lnx] = r2[1 : 1 + lmx][gx : gx + lnx]
        y[j0 : j0 + lny] = r2[1 + lmx : 1 + lmx + lmy][gy : gy + lny]
        z[k0 : k0 + lnz] = r2[1 + lmx + lmy : 1 + lmx + lmy + lmz][gz : gz + lnz]
    return A, t, x, y, z


def spectral_curl(A, d):
    """B = curl A on a periodic grid, computed with FFTs: exact for the Fourier interpolant of A,
    so no finite-difference error is added in post-processing.  The Nyquist wavenumber of an
    even-sized axis is set to zero (its derivative is undefined for a real field).

    Parameters: A: array (3, nx, ny, nz) indexed [component, x, y, z]
                d: grid spacings (dx, dy, dz)
    Returns:
                B: array (3, nx, ny, nz)
    """
    k = []
    for ax, (n, h) in enumerate(zip(A.shape[1:], d)):
        kk = 2 * np.pi * np.fft.fftfreq(n, d=h)
        if n % 2 == 0:
            kk[n // 2] = 0.0
        shape = [1, 1, 1]
        shape[ax] = n
        k.append(kk.reshape(shape))
    kx, ky, kz = k
    Ax, Ay, Az = np.fft.fftn(A, axes=(1, 2, 3))
    Bk = 1j * np.array([ky * Az - kz * Ay, kz * Ax - kx * Az, kx * Ay - ky * Ax])
    return np.real(np.fft.ifftn(Bk, axes=(1, 2, 3)))


def read_VAR_B(input_filename):
    """B = curl A from a Pencil snapshot, with the curl taken spectrally (spectral_curl).

    Parameters: input_filename: 'path/to/data/N.VAR' -> snapshot data/proc*/VARN
    Returns:
                Bx, By, Bz: arrays (nx, ny, nz) indexed [x, y, z];  Lbox: [Lx, Ly, Lz]
    """
    snapshotid = int(input_filename.split("/")[-1].split(".VAR")[0])
    datafolder = input_filename.split(str(snapshotid) + ".VAR")[0] or "."
    try:
        import pencil as pc
    except ImportError:
        pc = None
    if pc is not None:
        snap = pc.read.var(datadir=datafolder, ivar=snapshotid, trimall=True)
        aa = snap.aa.transpose(0, 3, 2, 1)  # (3, nz, ny, nx) -> (3, nx, ny, nz)
        t, x, y, z = snap.t, snap.x, snap.y, snap.z
        uu = snap.uu.transpose(3, 2, 1, 0)
        print("vrms", np.sqrt(np.mean(np.sum(uu * uu, axis=-1))))
    else:
        print("pencil package not found: reading the snapshot directly")
        aa, t, x, y, z = read_pencil_var(datafolder, f"VAR{snapshotid}")
    d = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
    Lbox = [len(x) * d[0], len(y) * d[1], len(z) * d[2]]
    print(f"Pencil VAR{snapshotid}: t = {float(t):.4f}, grid {aa.shape[1:]}, first point "
          f"x0 = {x[0]:.6f}, y0 = {y[0]:.6f} (pass -x0 to eigen_patterns.py if non-zero)")
    Bx, By, Bz = spectral_curl(aa, d)
    return Bx, By, Bz, Lbox

if __name__ == "__main__":

    import argparse as ap

    parser = ap.ArgumentParser(description="Generate ICs for ABC flow")

    parser.add_argument(
        "-res_z",
        "--grid_z_resolution",
        help="resolution of the grid (for SPH)",
        default=32,
        type=int,
    )
    parser.add_argument(
        "-res_xy",
        "--grid_xy_resolution",
        help="resolution of the grid (for SPH) in x and y",
        default=32,
        type=int,
    )
    parser.add_argument(
        "-in",
        "--input_filename",
        help="source file name",
        default="",
        type=str,
    )
    parser.add_argument(
        "-modez",
        help="index of the z Fourier mode, kz = 2*pi*modez/Lz (kz = 1 in Lz = 4*pi -> 2)",
        default=2,
        type=int,
    )
    parser.add_argument(
        "-out",
        "--out_slice",
        help="output slice name",
        default="slice.h5",
        type=str,
    )
    parser.add_argument(
        "-outim",
        "--out_image",
        help="output image name",
        default="im.png",
        type=str,
    )
    args = parser.parse_args()

    # Creating grid for B field and reading files
    res_z = args.grid_z_resolution
    res_xy = args.grid_xy_resolution

    filename = args.input_filename
    extension = filename.split(".")[-1]
    Bx, By, Bz = np.zeros([3, res_xy, res_xy, res_z])
    Lbox = [0, 0, 0]
    if extension == "hdf5":
        Bx, By, Bz, Lbox = read_hdf5_B(filename, res_xy=res_xy, res_z=res_z)
    #    elif extension == 'fits':
    #        Bx, By, Bz, Lbox = read_fits_B(filename)
    elif extension == 'VAR':
        # Pencil: the grid sets the resolution (-res_xy, -res_z are ignored)
        Bx, By, Bz, Lbox = read_VAR_B(filename)
        if Bx.shape[0] != Bx.shape[1]:
            raise Exception("Pencil grid must have nx == ny, got ", Bx.shape)
        res_xy, res_z = Bx.shape[0], Bx.shape[2]
    else:
        raise Exception(
            "input must be a SWIFT .hdf5 snapshot or a Pencil 'data/N.VAR', not ", filename
        )

    # Plotting

    # Pick specific mode, nz
    modez = args.modez

    # Load grid
    B_data = np.array([Bx, By, Bz])

    # Test plotter
    # B_data = generate_test_field(data)

    # Normalize B
    B2 = B_data[0][:, :, :] ** 2 + B_data[1][:, :, :] ** 2 + B_data[2][:, :, :] ** 2
    Brms = np.sqrt(np.mean(B2.flatten()))
    print("B_rms = ", Brms)
    B_data /= Brms

    # Extract modes
    k, b = transform_B_along_z(
        B_data[0], B_data[1], B_data[2], res_xy=res_xy, res_z=res_z
    )

    # Print mode momentum
    print("Mode: k=", k[modez])

    # Calculate |b_i|, <b_x>,  <b_y>
    b_abs = np.abs(b[:, :, modez, :])
    bx_mean = np.mean(b[:, :, modez, 0])
    # by_mean = np.mean(b[:,:,modez,1])
    bx_phase = bx_mean / np.abs(bx_mean)
    # by_phase = by_mean/np.abs(by_mean)

    # Re(<b_x>)>0
    # if np.real(bx_phase)<0:
    #    bx_phase*=-1
    # Re(<b_y>)>0
    # if np.real(by_phase)<0:
    #    by_phase*=-1

    # Display phase
    print("Complex phase arg(b_x):", np.angle(bx_phase))
    # print('Complex phase arg(b_y):',np.angle(by_phase))

    # No phase fixing
    theta_zero = np.angle(b[:, :, modez, :])

    # Im(<b_x>)=0
    theta = np.angle(b[:, :, modez, :] / bx_phase)
    # Im(<b_y>)=0
    # theta_tilde = np.angle(b[:,:,modez,:]/by_phase)

    # Fix phases to the interval \phi \in [-\pi,\pi]
    # theta_zero[theta_zero>np.pi] = np.pi - theta_zero[theta_zero>np.pi]
    # theta_zero[theta_zero<-np.pi] = theta_zero[theta_zero<-np.pi] - np.pi

    # theta[theta>np.pi] = np.pi - theta[theta>np.pi]
    # theta[theta<-np.pi] = theta[theta<-np.pi] - np.pi
    # theta_tilde[theta_tilde>np.pi] = np.pi - theta_tilde[theta_tilde>np.pi]
    # theta_tilde[theta_tilde<-np.pi] = theta_tilde[theta_tilde<-np.pi] - np.pi

    # Generate figure
    nx = 3
    ny = 2
    fig, ax = plt.subplots(ny, nx, sharey=True, sharex=True, figsize=(5 * nx, 5 * ny))

    from matplotlib.ticker import AutoMinorLocator

    for a in ax.flat:
        # Ticks on all four sides
        a.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            right=True,
            length=7,
            width=1.2,
            labelsize=16,
        )

        a.tick_params(
            axis="both",
            which="minor",
            direction="in",
            top=True,
            right=True,
            length=4,
            width=1.0,
        )

        # Add minor ticks
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())

    # Number of color levels
    nlev = 100  # 17

    # Colormaps
    # Custom colormaps
    # Select colormaps for r and theta
    cmapr = "plasma"  # c1r
    cmaptheta = "twilight_shifted"  # c1r #"RdBu_r"

    # Add plots
    a00 = ax[0, 0].pcolormesh(
        periodic_continue(b_abs[:, :, 0].T),
        cmap=cmapr,
        vmin=0,
        vmax=1.5,
        shading="nearest",  # important!
    )
    a01 = ax[0, 1].pcolormesh(
        periodic_continue(b_abs[:, :, 1].T),
        cmap=cmapr,
        vmin=0,
        vmax=1.5,
        shading="nearest",  # important!
    )
    a02 = ax[0, 2].pcolormesh(
        periodic_continue(b_abs[:, :, 2].T),
        cmap=cmapr,
        vmin=0,
        vmax=1.5,
        shading="nearest",  # important!
    )
    a10 = ax[1, 0].pcolormesh(
        periodic_continue(theta[:, :, 0].T),
        cmap=cmaptheta,
        vmin=-np.pi,
        vmax=np.pi,
        shading="nearest",  # important!
    )
    a11 = ax[1, 1].pcolormesh(
        periodic_continue(theta[:, :, 1].T),
        cmap=cmaptheta,
        vmin=-np.pi,
        vmax=np.pi,
        shading="nearest",  # important!
    )
    a12 = ax[1, 2].pcolormesh(
        periodic_continue(theta[:, :, 2].T),
        cmap=cmaptheta,
        vmin=-np.pi,
        vmax=np.pi,
        shading="nearest",  # important!
    )

    # Add ticks and labels
    locs = [res_xy / 4, res_xy / 2, 3 * res_xy / 4, res_xy]
    labels = [r"$\frac{\pi}{2}$", r"$\pi$", r"$\frac{3\pi}{2}$", r"$2\pi$"]
    # labels = [r'$\pi/2$', r'$\pi$', r'$3\pi/2$',r'$2\pi$']
    for ii in range(ny):
        # ax[ii, 0].set_ylabel(r"$y$",fontsize=20)
        ax[ii, 0].set_yticks(locs, labels, fontsize=25)
    for ii in range(ny):
        for jj in range(nx):
            # ax[ii, jj].set_xlabel(r"$x$",fontsize=20)
            ax[ii, jj].set_xticks(locs, labels, fontsize=25)
            # ax[ii, jj].set_aspect("equal")
            ax[ii, jj].set_aspect("equal", adjustable="box")
            ax[ii, jj].set_xlim(0, res_xy)
            ax[ii, jj].set_ylim(0, res_xy)

    ax[0, 0].tick_params(color="w", which="major")
    ax[0, 0].tick_params(color="w", which="minor")

    ax[0, 1].tick_params(color="w", which="major")
    ax[0, 1].tick_params(color="w", which="minor")

    ax[0, 2].tick_params(color="w", which="major")
    ax[0, 2].tick_params(color="w", which="minor")

    # Add colorbar
    ticks_r = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    ticks_phase = [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi]
    ticklabels_phase = [
        r"$-\pi$",
        r"$-\frac{\pi}{2}$",
        r"$0$",
        r"$\frac{\pi}{2}$",
        r"$\pi$",
    ]
    cbar1 = plt.colorbar(a00, ax=ax[0, 0], fraction=0.046, pad=0.04, ticks=ticks_r)
    cbar1.set_ticklabels(ticks_r, fontsize=15)
    # cbar1.set_label(r"$r_1 \;$",fontsize=20,rotation=0)
    ax[0, 0].set_title(r"$r_1 \;$", fontsize=30, pad=10)
    cbar2 = plt.colorbar(a01, ax=ax[0, 1], fraction=0.046, pad=0.04, ticks=ticks_r)
    cbar2.set_ticklabels(ticks_r, fontsize=15)
    # cbar2.set_label(r"$r_2 \;$",fontsize=20,rotation=0)
    ax[0, 1].set_title(r"$r_2 \;$", fontsize=30, pad=10)
    cbar3 = plt.colorbar(a02, ax=ax[0, 2], fraction=0.046, pad=0.04, ticks=ticks_r)
    cbar3.set_ticklabels(ticks_r, fontsize=15)
    # cbar3.set_label(r"$r_3 \;$",fontsize=20,rotation=0)
    ax[0, 2].set_title(r"$r_3 \;$", fontsize=30, pad=10)
    cbar4 = plt.colorbar(a10, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=ticks_phase)
    # cbar4.set_label(r"$\theta_1 \;$",fontsize=20, rotation=0)
    cbar4.set_ticklabels(ticklabels_phase, fontsize=20)
    ax[1, 0].set_title(r"$\theta_1 \;$", fontsize=30, pad=10)
    cbar5 = plt.colorbar(a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=ticks_phase)
    # cbar5.set_label(r"$\theta_2 \;$",fontsize=20,rotation=0)
    cbar5.set_ticklabels(ticklabels_phase, fontsize=20)
    ax[1, 1].set_title(r"$\theta_2 \;$", fontsize=30, pad=10)
    cbar6 = plt.colorbar(a12, ax=ax[1, 2], fraction=0.046, pad=0.04, ticks=ticks_phase)
    # cbar6.set_label(r"$\theta_3 \;$",fontsize=20,rotation=0)
    cbar6.set_ticklabels(ticklabels_phase, fontsize=20)
    ax[1, 2].set_title(r"$\theta_3 \;$", fontsize=30, pad=10)
    # cbar7 = plt.colorbar(a20, ax=ax[2, 0], fraction=0.046, pad=0.04, ticks=ticks_phase)
    # cbar7.set_label(r"$\tilde{\theta}_1(x,y) \;$")
    # cbar8 = plt.colorbar(a21, ax=ax[2, 1], fraction=0.046, pad=0.04, ticks=ticks_phase)
    # cbar8.set_label(r"$\tilde{\theta}_2(x,y) \;$")
    # cbar9 = plt.colorbar(a22, ax=ax[2, 2], fraction=0.046, pad=0.04, ticks=ticks_phase)
    # cbar9.set_label(r"$\tilde{\theta}_3(x,y) \;$")

    # save slices
    import h5py

    slice_dict = {
        "r1": periodic_continue(b_abs[:, :, 0].T),
        "r2": periodic_continue(b_abs[:, :, 1].T),
        "r3": periodic_continue(b_abs[:, :, 2].T),
        "theta1": periodic_continue(theta[:, :, 0].T),
        "theta2": periodic_continue(theta[:, :, 1].T),
        "theta3": periodic_continue(theta[:, :, 2].T),
        "ReB1": np.real(b[:, :, modez, 0] / bx_phase),
        "ImB1": np.imag(b[:, :, modez, 0] / bx_phase),
        "ReB2": np.real(b[:, :, modez, 1] / bx_phase),
        "ImB2": np.imag(b[:, :, modez, 1] / bx_phase),
        "ReB3": np.real(b[:, :, modez, 2] / bx_phase),
        "ImB3": np.imag(b[:, :, modez, 2] / bx_phase),
    }
    with h5py.File(args.out_slice, "w") as f:
        for key, arr in slice_dict.items():
            f.create_dataset(key, data=arr)

    fig.tight_layout()

    plt.savefig(args.out_image, dpi=220)
