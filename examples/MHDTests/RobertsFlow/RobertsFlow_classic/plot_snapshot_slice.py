from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import numpy as np
import sys
import h5py

filename = sys.argv[1]

slice_height = sys.argv[3]

projection = sys.argv[4]



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
    #dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    #dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]


data = load(filename)

# getting values from snapshots

img_res = 512

if projection=='yz':
    rotation_matrix = np.array([[0,1,0],[0,0,1],[1,0,0]])
    for i in range(len(data.gas.coordinates)):
        data.gas.coordinates[i] = np.matmul(rotation_matrix, data.gas.coordinates[i])
        data.gas.velocities[i] = np.matmul(rotation_matrix, data.gas.velocities[i])
        data.gas.magnetic_flux_densities[i] = np.matmul(rotation_matrix, data.gas.magnetic_flux_densities[i])
        data.gas.magnetic_vector_potentials[i] =np.matmul(rotation_matrix, data.gas.magnetic_vector_potentials[i])
if projection=='xz':
    rotation_matrix = np.array([[1,0,0],[0,0,1],[0,1,0]])
    for i in range(len(data.gas.coordinates)):
        data.gas.coordinates[i] = np.matmul(rotation_matrix, data.gas.coordinates[i])
        data.gas.velocities[i] = np.matmul(rotation_matrix, data.gas.velocities[i])
        data.gas.magnetic_flux_densities[i] = np.matmul(rotation_matrix, data.gas.magnetic_flux_densities[i])
        data.gas.magnetic_vector_potentials[i] =np.matmul(rotation_matrix, data.gas.magnetic_vector_potentials[i])


x = data.gas.coordinates[:, 0].value
y = data.gas.coordinates[:, 1].value
z = data.gas.coordinates[:, 2].value
rho = data.gas.densities.value
h = data.gas.smoothing_lengths.value
v = data.gas.velocities.value
P = data.gas.pressures.value
B = data.gas.magnetic_flux_densities.value
A = data.gas.magnetic_vector_potentials.value
divB = data.gas.magnetic_divergences.value

zmin = np.min(z)
zmax = np.max(z)

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
data.gas.mass_weighted_Ax = data.gas.masses * A[:, 0]
data.gas.mass_weighted_Ay = data.gas.masses * A[:, 1]
data.gas.mass_weighted_Az = data.gas.masses * A[:, 2]
data.gas.mass_weighted_divB = data.gas.masses * divB

def make_slice(key, slice_frac_z=float(slice_height)):
        res = slice_gas(
        data,
        z_slice = slice_frac_z *data.metadata.boxsize[2],
        resolution=img_res,
        project=key,
        parallel=True,
        periodic=True,
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
A = np.zeros((l, 3))

rho = make_slice("mass_weighted_rho").value.flatten() / mass_map
h = make_slice("mass_weighted_h").value.flatten() / mass_map
v[:, 0] = make_slice("mass_weighted_vx").value.flatten() / mass_map
v[:, 1] = make_slice("mass_weighted_vy").value.flatten() / mass_map
v[:, 2] = make_slice("mass_weighted_vz").value.flatten() / mass_map
P = make_slice("mass_weighted_P").value.flatten() / mass_map
B[:, 0] = make_slice("mass_weighted_Bx").value.flatten() / mass_map
B[:, 1] = make_slice("mass_weighted_By").value.flatten() / mass_map
B[:, 2] = make_slice("mass_weighted_Bz").value.flatten() / mass_map

A[:, 0] = make_slice("mass_weighted_Ax").value.flatten() / mass_map
A[:, 1] = make_slice("mass_weighted_Ay").value.flatten() / mass_map
A[:, 2] = make_slice("mass_weighted_Az").value.flatten() / mass_map

divB = make_slice("mass_weighted_divB").value.flatten() / mass_map



bb = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
aa = np.sqrt(A[:, 0] ** 2 + A[:, 1] ** 2 + A[:, 2] ** 2)
Pmag = (B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2) / 2
vv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)

R0 = np.abs(divB) * h / (bb + np.abs(divB) * h)

reg_err = 1.001 * 1e-2
err_upper_bnd = 0.9
above_noise = 10

mask_lower_bound = R0<reg_err
R0[mask_lower_bound]=reg_err
mask_upper_bound = R0>err_upper_bnd
R0[mask_upper_bound]=err_upper_bnd


from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import ticker


# plot everything
fig, ax = plt.subplots(1,3, sharex=True, figsize=(18, 5))

#if project == "xy":
new_x = np.linspace(np.min(x), np.max(x), dimx)
new_y = np.linspace(np.min(y), np.max(y), dimy)
#elif project == "yz":
#    new_x = np.linspace(np.min(y.value), np.max(y.value), dimx)
#    new_y = np.linspace(np.min(z.value), np.max(z.value), dimy)
#elif project == "zx":
#    new_x = np.linspace(np.min(z.value), np.max(z.value), dimx)
#    new_y = np.linspace(np.min(x.value), np.max(x.value), dimy)


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
        to_plot = ax[j].contourf(
            new_x,
            new_y,
            Q.transpose(),
            levels=np.array(levels),
            locator=ticker.LogLocator(),
            cmap=cmap,
        )
    else:
        to_plot = ax[j].contourf(
            new_x, new_y, Q.transpose(), levels=np.array(levels), cmap=cmap
        )

    fig.colorbar(to_plot, ticks=levels_short)
    ax[j].set_ylabel(Q_name)
    ax[j].set_xlim(min(new_x),max(new_x))
    ax[j].set_ylim(min(new_y),max(new_y))
    return 0

def make_slice_plot(Q, cmin, cmax, i, j, Q_name):
    slice_Q = Q[:, int(len(Q) / 2)]
    ax[j].plot(x, slice_Q)
    ax[j].plot(x, max(slice_Q) * np.ones(len(slice_Q)), "--")
    ax[j].set_ylim([cmin, cmax])
    ax[j].set_yscale("log")
    ax[j].set_ylabel(Q_name)
    return 0


Pmag = bb

make_density_plot(
    bb.reshape((dimx, dimy)),
    np.min(bb),
    1.01*np.max(bb),
    0,
    0,
    "B",
    c_res=100,
    log_sc=False
)
ax[0].streamplot(new_x, new_y, np.transpose(B[:,0].reshape((dimx, dimy))), np.transpose(B[:,1].reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)

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
ax[1].streamplot(new_x, new_y, np.transpose(v[:,0].reshape((dimx, dimy))),np.transpose(v[:,1].reshape((dimx, dimy))), color='w', density=1.0, linewidth=0.5, arrowsize=0.8)

#make_density_plot(
#    R0.reshape((dimx, dimy)),
#    np.min(R0),
#    1.01*np.max(R0),
#    0,
#    2,
#    "R0",
#    c_res=100,
#    log_sc=False
#)
#ax[2].streamplot(new_x, new_y, np.transpose(B[:,0].reshape((dimx, dimy))), np.transpose(B[:,1].reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)

make_density_plot(
    aa.reshape((dimx, dimy)),
    np.min(aa),
    1.01*np.max(aa),
    0,
    2,
    "|A|",
    c_res=100,
    log_sc=False
)
ax[2].streamplot(new_x, new_y, np.transpose(A[:,0].reshape((dimx, dimy))), np.transpose(A[:,1].reshape((dimx, dimy))), color='w', density=2.0, linewidth=0.5, arrowsize=0.8)

ax[1].set_title(f"t={t:.2e}")
ax[0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
fig.tight_layout()

plt.savefig(sys.argv[2], dpi=70)
