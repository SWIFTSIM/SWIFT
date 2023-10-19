from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import numpy as np
import sys
import h5py

filename = sys.argv[1]

with h5py.File(filename, "r") as handle:
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"][0]
    git = handle["Code"].attrs["Git Revision"]
    gitBranch = handle["Code"].attrs["Git Branch"]
    scheme = handle["/HydroScheme"].attrs["Scheme"]
    kernel = handle["/HydroScheme"].attrs["Kernel function"]
    neighbours = handle["/HydroScheme"].attrs["Kernel target N_ngb"]
    mu0 = handle["/PhysicalConstants/InternalUnits"].attrs["vacuum_permeability"]
    dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]


filename = sys.argv[1]
data = load(filename)

# getting values from snapshots

img_res = 1000
project = "xy"

x = data.gas.coordinates[:, 0]
y = data.gas.coordinates[:, 1]
z = data.gas.coordinates[:, 2]
rho = data.gas.densities.value
h = data.gas.smoothing_lengths.value
v = data.gas.velocities.value
P = data.gas.pressures.value
B = data.gas.magnetic_flux_densities.value

data.gas.mass_weighted_x = data.gas.masses * x[:]
data.gas.mass_weighted_y = data.gas.masses * y[:]
data.gas.mass_weighted_z = data.gas.masses * z[:]
data.gas.mass_weighted_rho = data.gas.masses * rho[:]
data.gas.mass_weighted_h = data.gas.masses * h[:]
data.gas.mass_weighted_vx = data.gas.masses * v[:, 0]
data.gas.mass_weighted_vy = data.gas.masses * v[:, 1]
data.gas.mass_weighted_vz = data.gas.masses * v[:, 2]
data.gas.mass_weighted_P = data.gas.masses * P[:]
data.gas.mass_weighted_Bx = data.gas.masses * B[:, 0]
data.gas.mass_weighted_By = data.gas.masses * B[:, 1]
data.gas.mass_weighted_Bz = data.gas.masses * B[:, 2]


def make_slice(key):
    res = slice_gas(
        data,
        z_slice=0.01 * data.metadata.boxsize[2],
        resolution=img_res,
        project=key,
        parallel=True,
    )
    return res


divreg = 1e-30

masses = make_slice("masses").value
dimy = len(masses)
dimx = len(masses[0])
mass_map = masses.flatten()
mass_map = mass_map + divreg

l = len(mass_map)

v = np.zeros((l, 3))
B = np.zeros((l, 3))

rho = make_slice("mass_weighted_rho").value.flatten() / mass_map
h = make_slice("mass_weighted_h").value.flatten() / mass_map
v[:, 0] = make_slice("mass_weighted_vx").value.flatten() / mass_map
v[:, 1] = make_slice("mass_weighted_vy").value.flatten() / mass_map
v[:, 2] = make_slice("mass_weighted_vz").value.flatten() / mass_map
P = make_slice("mass_weighted_P").value.flatten() / mass_map
B[:, 0] = make_slice("mass_weighted_Bx").value.flatten() / mass_map
B[:, 1] = make_slice("mass_weighted_By").value.flatten() / mass_map
B[:, 2] = make_slice("mass_weighted_Bz").value.flatten() / mass_map

bb = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
Pmag = (B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2) / 2
vv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)

reg_err = 1.001 * 1e-2
err_upper_bnd = 0.9
above_noise = 10

from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import ticker


# plot everything
fig, ax = plt.subplots(2, 2, sharex=True, figsize=(14, 10))

if project == "xy":
    new_x = np.linspace(np.min(x.value), np.max(x.value), dimx)
    new_y = np.linspace(np.min(y.value), np.max(y.value), dimy)
elif project == "yz":
    new_x = np.linspace(np.min(y.value), np.max(y.value), dimx)
    new_y = np.linspace(np.min(z.value), np.max(z.value), dimy)
elif project == "zx":
    new_x = np.linspace(np.min(z.value), np.max(z.value), dimx)
    new_y = np.linspace(np.min(x.value), np.max(x.value), dimy)


def make_color_levels(cmin, cmax, c_res=10, log_sc=True):
    cmax += divreg
    if log_sc:
        levmin = int(np.floor(np.log10(cmin)))
        levmax = int(np.ceil(np.log10(cmax)))
        levels = []
        levels_short = []
        for i in range(levmin, levmax):
            levels_short += [10 ** i]
            for j in range(int(c_res / 10), c_res):
                levels += [(10 / c_res * j) * 10 ** i]

    else:
        levels = [cmin + (cmax - cmin) / c_res * k for k in range(c_res)]
        levels_short = [cmin + (cmax - cmin) / c_res * k for k in range(0, c_res, 10)]
    return levels, levels_short


def make_density_plot(
    Q, cmin, cmax, i, j, Q_name, c_res=10, log_sc=True, cmap="viridis"
):
    levels, levels_short = make_color_levels(cmin, cmax, c_res, log_sc)
    if log_sc:
        to_plot = ax[i][j].contourf(
            new_x,
            new_y,
            Q.transpose(),
            levels=np.array(levels),
            locator=ticker.LogLocator(),
            cmap=cmap,
        )
    else:
        to_plot = ax[i][j].contourf(
            new_x, new_y, Q.transpose(), levels=np.array(levels), cmap=cmap
        )
    fig.colorbar(to_plot, ticks=levels_short)
    ax[i][j].set_ylabel(Q_name)
    return 0


def make_slice_plot(Q, cmin, cmax, i, j, Q_name):
    slice_Q = Q[:, int(len(Q) / 2)]
    ax[i][j].plot(x, slice_Q)
    ax[i][j].plot(x, max(slice_Q) * np.ones(len(slice_Q)), "--")
    ax[i][j].set_ylim([cmin, cmax])
    ax[i][j].set_yscale("log")
    ax[i][j].set_ylabel(Q_name)
    return 0


Pmag = bb

make_density_plot(
    Pmag.reshape((dimx, dimy)),
    np.min(Pmag),
    1.1 * np.max(Pmag),
    0,
    0,
    "B",
    c_res=100,
    log_sc=False,
)
make_density_plot(
    P.reshape((dimx, dimy)),
    np.min(P),
    1.1 * np.max(P),
    1,
    0,
    "Pressure",
    c_res=100,
    log_sc=False,
)
make_density_plot(
    vv.reshape((dimx, dimy)),
    np.min(vv),
    1.1 * np.max(vv),
    0,
    1,
    "Velocity",
    c_res=100,
    log_sc=False,
)
make_density_plot(
    rho.reshape((dimx, dimy)),
    np.min(rho),
    1.1 * np.max(rho),
    1,
    1,
    "Density",
    c_res=100,
    log_sc=False,
)

ax[0, 1].set_title(f"t={t:.2e}")
ax[0, 0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
fig.tight_layout()

plt.savefig(sys.argv[2], dpi=300)
