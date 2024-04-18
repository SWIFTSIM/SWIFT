from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
import numpy as np
import sys
import h5py
from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter
from matplotlib import colors

filename = sys.argv[1]

slice_height = sys.argv[3]

projection = sys.argv[4]

prefered_color = "magma"

cpts = 100


to_plot = "|B|_rho_|v|_R0"  # 'B' or 'A' or 'errors' or |B|_rho_|v|_R0

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
    # dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    # dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]


data = load(filename)

# getting values from snapshots

img_res = 512
divreg = 1e-30
reg_err = 0.01

if projection == "yz":
    rotation_matrix = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    for i in range(len(data.gas.coordinates)):
        data.gas.coordinates[i] = np.matmul(rotation_matrix, data.gas.coordinates[i])
        data.gas.velocities[i] = np.matmul(rotation_matrix, data.gas.velocities[i])
        data.gas.magnetic_flux_densities[i] = np.matmul(
            rotation_matrix, data.gas.magnetic_flux_densities[i]
        )
        data.gas.magnetic_vector_potentials[i] = np.matmul(
            rotation_matrix, data.gas.magnetic_vector_potentials[i]
        )
if projection == "xz":
    rotation_matrix = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    for i in range(len(data.gas.coordinates)):
        data.gas.coordinates[i] = np.matmul(rotation_matrix, data.gas.coordinates[i])
        data.gas.velocities[i] = np.matmul(rotation_matrix, data.gas.velocities[i])
        data.gas.magnetic_flux_densities[i] = np.matmul(
            rotation_matrix, data.gas.magnetic_flux_densities[i]
        )
        data.gas.magnetic_vector_potentials[i] = np.matmul(
            rotation_matrix, data.gas.magnetic_vector_potentials[i]
        )


def abs_vec(vec):
    res = np.sqrt(vec[:, 0] ** 2 + vec[:, 1] ** 2 + vec[:, 2] ** 2)
    return res


def dot_vec(vec1, vec2):
    res = vec1[:, 0] * vec2[:, 0] + vec1[:, 1] * vec2[:, 1] + vec1[:, 2] * vec2[:, 2]
    return res


def cross_vec(vec1, vec2):
    res_vec = np.zeros((len(vec1), 3))
    res_vec[:, 0] = vec1[:, 1] * vec2[:, 2] - vec1[:, 2] * vec2[:, 1]
    res_vec[:, 1] = vec1[:, 2] * vec2[:, 0] - vec1[:, 0] * vec2[:, 2]
    res_vec[:, 2] = vec1[:, 0] * vec2[:, 1] - vec1[:, 1] * vec2[:, 0]
    return res_vec


def rms_vec(vec):
    res = np.sqrt(
        np.mean(vec[:, 0] * vec[:, 0] + vec[:, 1] * vec[:, 1] + vec[:, 2] * vec[:, 2])
    )
    return res


# see available fields in snapshot
# print(data.metadata.gas_properties.field_names)
print(data.metadata)
a = data.metadata.cosmology["Scale-factor"][0]
z_to_display = np.round(data.metadata.cosmology["Redshift"][0], 1)
box_size_physical = a * data.metadata.boxsize
# print('Scale factor, a=',a)
print("z = ", z_to_display)

# Get physical quantities
x = data.gas.coordinates[:, 0].to_physical()
y = data.gas.coordinates[:, 1].to_physical()
z = data.gas.coordinates[:, 2].to_physical()
# print(x.cosmo_factor.expr)
# print(x)
print(x[0])
# print(x.to_physical())

rho = data.gas.densities
h = data.gas.smoothing_lengths
v = data.gas.velocities
P = data.gas.pressures
B = data.gas.magnetic_flux_densities

R0 = data.gas.r0
R1 = data.gas.r1
R2 = data.gas.r2
R3 = data.gas.r3

print(x[0], B[0], rho[0])

print("Average error metrics:")

print("R0 = ", np.log10(np.mean(R0.to_physical())))
print("R1 = ", np.log10(np.mean(R1.to_physical())))
print("R2 = ", np.log10(np.mean(R2.to_physical())))
print("R3 = ", np.log10(np.mean(R3.to_physical())))


# Get RMS values
Brms = rms_vec(B).to_physical()
vrms = rms_vec(v).to_physical()

if to_plot == "A":
    A = data.gas.magnetic_vector_potentials
    Arms = rms_vec(A.value).to_physical()

# Get mean density
rho_mean = np.sum(rho[:]).to_physical() / len(rho)

# Scan for particles with maximal magnetic fields and display
Babs = abs_vec(B).to_physical()
indx = np.argmax(Babs)
z_slice_frac = np.round(z[indx].value / box_size_physical.value[2], 4)
Babs_particle = Babs[indx] / Brms
h_dev_Lbox = h.to_physical()[indx] / box_size_physical

print(
    "Maximal MF for particle #"
    + str(indx)
    + " at z_slice/Lbox = "
    + str(z_slice_frac)
    + " with |B|/Brms="
    + str(Babs_particle.value)
    + " and h/L_box="
    + str(h_dev_Lbox)
    + " at x="
    + str(x[indx])
    + ", y="
    + str(y[indx])
)

# Scan for place with maximal error and display
maxR0 = np.max(R0.value)
if maxR0 >= 1:
    indx = np.argmax(R0)
    z_slice_frac = np.round(z[indx].value / box_size_physical.value[2], 4)

    Babs_particle = Babs[indx] / Brms
    print(
        "Warning: there are particles with lagre divB error. Maximal error for particle #"
        + str(np.argmax(R0))
        + " at z_slice/Lbox = "
        + str(z_slice_frac)
        + " with |B|/Brms="
        + str(Babs_particle.value)
    )
else:
    print("R0 errors are acceptable")


# function to prepare colorbar tick values
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


########################################## Make histogram
if to_plot == "hist":
    min_metric = 0.1
    mask = R0.value >= min_metric
    quantity = (
        abs_vec(B).to_physical() / Brms
    ).value  # (rho.to_physical()/rho_mean).value#(abs_vec(v).to_physical()/vrms).value
    qname = "$|B|/B_{rms}$"  #'rho/<rho>' #'$|v|/v_{rms}$'
    Nbins = 100
    # bins = np.linspace(np.min(quantity),np.max(quantity),Nbins)
    bins = np.logspace(
        int(np.log10(np.min(quantity))) - 1, int(np.log10(np.max(quantity))) + 1, Nbins
    )
    fig, ax = plt.subplots()
    plt.hist(quantity, bins=bins, color="blue", label=qname, alpha=0.5, density=False)
    plt.hist(
        quantity[mask],
        bins=bins,
        color="red",
        label=qname + "($R_0$>" + str(min_metric) + ")",
        alpha=0.5,
        density=False,
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    ax.set_ylabel("$N_{particles}$/bin")
    ax.set_xlabel(qname)
    ax.set_title(f"z={z_to_display:}")
    plt.grid()
    plt.savefig(sys.argv[2], dpi=100)
    exit()
########################################## Make correlation histograms


def prepare_histogram_density_plot(data_x, data_y, Nbins):
    bins_x = np.logspace(
        int(np.log10(np.min(data_x))) - 1, int(np.log10(np.max(data_x))) + 1, Nbins
    )
    bins_y = np.logspace(
        int(np.log10(np.min(data_y))) - 1, int(np.log10(np.max(data_y))) + 1, Nbins
    )
    hist, xedges, yedges = np.histogram2d(data_x, data_y, bins=[bins_x, bins_y])
    # building uncorrelated histogram
    hist_x, _ = np.histogram(data_x, bins=bins_x)
    hist_y, _ = np.histogram(data_y, bins=bins_y)
    hist_uncorr = np.outer(hist_x, hist_y).astype("float64")

    # norm
    min_level = 1 / np.sum(hist)
    hist /= np.sum(hist)
    hist_uncorr /= np.sum(hist_uncorr)

    # print(np.sum(hist),np.sum(hist_uncorr))
    # regularize
    hist[hist < min_level] = min_level
    hist_uncorr[hist_uncorr < min_level] = min_level

    # build correlation_funciton
    correlation_funcion = hist - hist_uncorr

    X, Y = np.meshgrid(xedges, yedges)
    return X, Y, hist, hist_uncorr, correlation_funcion


def plot_histogram_density_plot(quantity_x, quantity_y, qname_x, qname_y, Nbins=100):
    quantity_x = quantity_x.value
    quantity_y = quantity_y.value

    X, Y, hist, hist_uncorr, correlation_funcion = prepare_histogram_density_plot(
        quantity_x, quantity_y, Nbins
    )

    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    fig.suptitle(f"({qname_x},{qname_y}) space at z={z_to_display:}")

    ax[0].set_title(r"$f_{12}(x,y)$")
    ax[0].set_ylabel(qname_y)
    ax[0].set_xlabel(qname_x)
    ax[0].set_box_aspect(1)

    pax0 = ax[0].pcolormesh(X, Y, hist, cmap=prefered_color, norm=LogNorm())
    ax[0].grid()
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    fig.colorbar(pax0)

    ax[1].set_box_aspect(1)
    ax[1].set_title(r"$f_1(x)\cdot f_2(y)$")
    pax1 = ax[1].pcolormesh(X, Y, hist_uncorr, cmap=prefered_color, norm=LogNorm())
    ax[1].grid()
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    fig.colorbar(pax1)

    vmax = np.max(np.max(hist))
    corr_significance = np.sum(np.abs(correlation_funcion))
    divnorm = colors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    ax[2].set_box_aspect(1)
    ax[2].set_title(r"$G_{corr}(x,y)$")
    pax2 = ax[2].pcolormesh(X, Y, correlation_funcion, cmap="seismic", norm=divnorm)
    ax[2].set_xscale("log")
    ax[2].set_yscale("log")
    ax[2].grid()
    fig.colorbar(pax2)

    fig.tight_layout()

    plt.savefig(sys.argv[2], dpi=100)


if to_plot == "correlation_hist":
    q_x = rho.to_physical() / rho_mean
    q_y = R0  # abs_vec(B).to_physical()/Brms
    mask = q_y.value >= 0.0001
    plot_histogram_density_plot(
        q_x[mask], q_y[mask], r"$\rho/\langle \rho \rangle$", "$R_0$", Nbins=100
    )
    exit()


#################################################### Make slices


def make_slice(key, slice_frac_z=float(slice_height)):
    res = slice_gas(
        data,
        z_slice=slice_frac_z * data.metadata.boxsize[2],
        resolution=img_res,
        project=key,
        parallel=True,
        periodic=True,
    )
    return res


# Slicing preparation
masses = make_slice("masses")
mass_map = masses.flatten()
mass_map = mass_map + divreg * masses.units

l = len(mass_map)


def prepare_sliced_quantity(quantity, isvec=False):
    mass = data.gas.masses
    if isvec:
        # Prepare vector quantity for slicing
        data.gas.mass_weighted_temp_qx = mass * quantity[:, 0]
        data.gas.mass_weighted_temp_qy = mass * quantity[:, 1]
        data.gas.mass_weighted_temp_qz = mass * quantity[:, 2]

        # Apply slicing for each component
        sliced_quantity_x = make_slice("mass_weighted_temp_qx").flatten() / mass_map
        sliced_quantity_y = make_slice("mass_weighted_temp_qy").flatten() / mass_map
        sliced_quantity_z = make_slice("mass_weighted_temp_qz").flatten() / mass_map

        # Convert vector quantity to physical
        sliced_quantity_x = sliced_quantity_x.to_physical()
        sliced_quantity_y = sliced_quantity_y.to_physical()
        sliced_quantity_z = sliced_quantity_z.to_physical()

        # Join components together (should be optimized!!!)
        the_units = sliced_quantity_x.units
        sliced_quantity = np.stack(
            (sliced_quantity_x.value, sliced_quantity_y.value, sliced_quantity_z.value),
            axis=-1,
        )
        sliced_quantity = sliced_quantity * the_units
    else:
        # prepare scalar quantity for slicing
        data.gas.mass_weighted_temp_q = data.gas.masses * quantity[:]
        # apply slicing
        sliced_quantity = make_slice("mass_weighted_temp_q").flatten() / mass_map
        # Convert scalar quantity to physical
        sliced_quantity = sliced_quantity.to_physical()
    return sliced_quantity


# slicing scalars
rho = prepare_sliced_quantity(rho)
h = prepare_sliced_quantity(h)
P = prepare_sliced_quantity(P)

# slicing vectors
v = prepare_sliced_quantity(v, isvec=True)
B = prepare_sliced_quantity(B, isvec=True)
if to_plot == "A":
    A = prepare_sliced_quantity(A, isvec=True)

# Get sliced error metrics
R0 = prepare_sliced_quantity(R0)
R1 = prepare_sliced_quantity(R1)
R2 = prepare_sliced_quantity(R2)
R3 = prepare_sliced_quantity(R3)

# mask error metrics
reg_err = reg_err * R0.units
R0[R0 < reg_err] = reg_err
R1[R1 < reg_err] = reg_err
R2[R2 < reg_err] = reg_err
R3[R3 < reg_err] = reg_err

################################## plot slices

# Create XY grid for plotting

l = len(mass_map)
dimy = len(masses)
dimx = len(masses[0])
new_x = np.linspace(0.0, box_size_physical.value[0], dimx)
new_y = np.linspace(0.0, box_size_physical.value[1], dimy)


from matplotlib.ticker import FormatStrFormatter

# function to make single density plot of some quantity Q in range from cmin to cmax
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
            extend="both",
        )
    else:
        to_plot = ax[j].contourf(
            new_x,
            new_y,
            Q.transpose(),
            levels=np.array(levels),
            cmap=cmap,
            extend="both",
        )

    fig.colorbar(to_plot, ticks=levels_short)
    ax[j].set_ylabel(Q_name)
    ax[j].set_xlim(min(new_x), max(new_x))
    ax[j].set_ylim(min(new_y), max(new_y))
    ax[j].yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax[j].xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    return 0


# Plot vector potential
if to_plot == "A":

    # get quantities to plot

    Ax = A[:, 0] / Arms
    Ay = A[:, 1] / Arms
    Az = A[:, 2] / Arms

    # save units

    # Convert quantities to numbers before plotting
    Ax = Ax.value
    Ay = Ay.value
    Az = Az.value

    fig, ax = plt.subplots(1, 3, sharex=True, figsize=(6 * 3, 5))
    make_density_plot(
        Ax.reshape((dimx, dimy)),
        -1.0,
        1.0,
        0,
        0,
        "$A_x$/$A_{rms}$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        Ay.reshape((dimx, dimy)),
        -1.0,
        1.0,
        0,
        1,
        "$A_y$/$A_{rms}$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        Az.reshape((dimx, dimy)),
        -1.0,
        1.0,
        0,
        2,
        "$A_z$/$A_{rms}$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    ax[0].streamplot(
        new_x,
        new_y,
        np.transpose(Ax.reshape((dimx, dimy))),
        np.transpose(Ay.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.5,
        arrowsize=0.8,
    )
    ax[1].streamplot(
        new_x,
        new_y,
        np.transpose(Ax.reshape((dimx, dimy))),
        np.transpose(Ay.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.5,
        arrowsize=0.8,
    )
    ax[2].streamplot(
        new_x,
        new_y,
        np.transpose(Ax.reshape((dimx, dimy))),
        np.transpose(Ay.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.5,
        arrowsize=0.8,
    )
    ax[0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
    ax[1].set_title(f"z={z_to_display:}")
    ax[2].set_title("$z_{slice}/L_{box}$=" + slice_height)
    fig.tight_layout()
    plt.savefig(sys.argv[2], dpi=100)


# Plot magnetic fields
if to_plot == "B":

    # get quantities to plot

    Bx = B[:, 0] / Brms
    By = B[:, 1] / Brms
    Bz = B[:, 2] / Brms

    # save units

    # Convert quantities to numbers before plotting
    Bx = Bx.value
    By = By.value
    Bz = Bz.value

    fig, ax = plt.subplots(1, 3, sharex=True, figsize=(6 * 3, 5))
    make_density_plot(
        Bx.reshape((dimx, dimy)),
        -1.0,
        1.0,
        0,
        0,
        "$B_x$/$B_{rms}$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        By.reshape((dimx, dimy)),
        -1.0,
        1.0,
        0,
        1,
        "$B_y$/$B_{rms}$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        Bz.reshape((dimx, dimy)),
        -1.0,
        1.0,
        0,
        2,
        "$B_z$/$B_{rms}$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    ax[0].streamplot(
        new_x,
        new_y,
        np.transpose(Bx.reshape((dimx, dimy))),
        np.transpose(By.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.5,
        arrowsize=0.8,
    )
    ax[1].streamplot(
        new_x,
        new_y,
        np.transpose(Bx.reshape((dimx, dimy))),
        np.transpose(By.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.5,
        arrowsize=0.8,
    )
    ax[2].streamplot(
        new_x,
        new_y,
        np.transpose(Bx.reshape((dimx, dimy))),
        np.transpose(By.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.5,
        arrowsize=0.8,
    )
    ax[0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
    ax[1].set_title(f"z={z_to_display:}")
    ax[2].set_title("$z_{slice}/L_{box}$=" + slice_height)
    fig.tight_layout()
    plt.savefig(sys.argv[2], dpi=100)


# Plot combination of Absolute values of B, v, plot density and classical error metric
if to_plot == "|B|_rho_|v|_R0":

    # get quantities to plot
    Babs = abs_vec(B) / Brms
    vabs = abs_vec(v) / vrms
    rho = rho / rho_mean
    Bx = B[:, 0]
    By = B[:, 1]
    vx = v[:, 0]
    vy = v[:, 1]

    # save units
    B_units = Babs.units
    v_units = vabs.units
    rho_units = rho.units

    # convert quantities to numbers before plotting
    Babs = Babs.value
    rho = rho.value
    vabs = vabs.value
    R0 = R0.value
    Bx = Bx.value
    By = By.value
    vx = vx.value
    vy = vy.value

    fig, ax = plt.subplots(1, 4, sharex=True, figsize=(6 * 4, 5))

    make_density_plot(
        Babs.reshape((dimx, dimy)),
        1e-4,
        1e2,
        0,
        0,
        "$|B|/B_{rms}$",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        rho.reshape((dimx, dimy)),
        1e-4,
        1e2,
        0,
        1,
        " $rho$ / <rho>",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        vabs.reshape((dimx, dimy)),
        1e-1,
        1e1,
        0,
        2,
        "$|v|/v_{rms}$",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        R0.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        3,
        "$R_0$",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    ax[0].streamplot(
        new_x,
        new_y,
        np.transpose(Bx.reshape((dimx, dimy))),
        np.transpose(By.reshape((dimx, dimy))),
        color="w",
        density=4.0,
        linewidth=0.1,
        arrowsize=0.2,
    )
    ax[2].streamplot(
        new_x,
        new_y,
        np.transpose(vx.reshape((dimx, dimy))),
        np.transpose(vy.reshape((dimx, dimy))),
        color="w",
        density=2.0,
        linewidth=0.25,
        arrowsize=0.4,
    )
    ax[0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
    ax[1].set_title(f"z={z_to_display:}")
    ax[2].set_title("$z_{slice}/L_{box}$=" + slice_height)
    fig.tight_layout()
    plt.savefig(sys.argv[2], dpi=300)


# Plot all error metrics
if to_plot == "errors":

    # save units

    # convert everything to numbers before plotting
    R0 = R0.value
    R1 = R1.value
    R2 = R2.value
    R3 = R3.value

    fig, ax = plt.subplots(
        1, 4, sharex=True, figsize=(24, 5)
    )  # for 3 plts 18 for 4 plts use 24

    make_density_plot(
        R0.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        0,
        "$R_0$",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        R1.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        1,
        "$R_1$",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        R2.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        2,
        "$R_2$",
        c_res=cpts,
        # log_sc=False,
        cmap=prefered_color,
    )
    make_density_plot(
        R3.reshape((dimx, dimy)),
        reg_err,
        1.0,
        0,
        3,
        "$R_3$",
        c_res=cpts,
        log_sc=False,
        cmap=prefered_color,
    )
    ax[0].set_title(f"Nneigh={int(neighbours[0]):}, Npart={len(data.gas.coordinates):}")
    ax[1].set_title(f"z={z_to_display:}")
    ax[2].set_title("$z_{slice}/L_{box}$=" + slice_height)
    fig.tight_layout()
    plt.savefig(sys.argv[2], dpi=100)
