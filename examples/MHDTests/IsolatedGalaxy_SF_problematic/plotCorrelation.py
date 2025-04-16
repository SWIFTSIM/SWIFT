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
import unyt

filename = sys.argv[1]

# slice_height = sys.argv[3]

# projection = sys.argv[4]

prefered_color = "magma"

cpts = 100


to_plot = "correlation_hist"  # 'B' or 'A' or 'errors' or |B|_rho_|v|_R0

data = load(filename)

time = data.metadata.time
time.convert_to_units(unyt.Gyr)
print("Plotting at %2f Gyr" % time)

# getting values from snapshots

img_res = 512
divreg = 1e-30
reg_err = 0.01

# see available fields in snapshot
# print(data.metadata.gas_properties.field_names)
# Get physical quantities
x = data.gas.coordinates[:, 0].to_physical()
y = data.gas.coordinates[:, 1].to_physical()
z = data.gas.coordinates[:, 2].to_physical()

rho = data.gas.densities
nH = rho.to(unyt.g / unyt.cm ** 3) / (1.67e-24 * unyt.g)
h = data.gas.smoothing_lengths
v = data.gas.velocities
P = data.gas.pressures
T = data.gas.temperatures


# Retrieve some information about the simulation run
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]


########################################## Make histogram
if to_plot == "hist":
    # min_metric = 0.1
    quantity = (
        nH  # T.to(unyt.K)#rho.to(unyt.g/unyt.cm**3)
    ).value  # (rho.to_physical()/rho_mean).value#(abs_vec(v).to_physical()/vrms).value
    qname = r"$n_H$ $[cm^{-3}]$"  # r"$\rho$ $[g/cm^3]$"  #'rho/<rho>' #'$|v|/v_{rms}$'
    Nbins = 100
    # bins = np.linspace(np.min(quantity),np.max(quantity),Nbins)
    bins = np.logspace(
        int(np.log10(np.min(quantity))) - 1, int(np.log10(np.max(quantity))) + 1, Nbins
    )
    fig, ax = plt.subplots()

    nx = 2
    ny = 1
    fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5 * nx, 5 * ny))

    ax[0].hist(quantity, bins=bins, color="blue", label=qname, alpha=0.5, density=False)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].legend()
    ax[0].set_ylabel("$N_{particles}$/bin")
    ax[0].set_xlabel(qname)
    ax[0].set_ylim([1, 1e5])
    # ax.set_title(f"z={z_to_display:}")
    ax[0].grid()

    # add panel with infromation about the run
    Np = len(quantity)
    text_common_args = dict(
        fontsize=10, ha="center", va="center", transform=ax[1].transAxes
    )

    ax[1].text(
        0.5,
        0.8,
        "Cooling halo with spin at time $t=%.2f$ Gyr" % data.metadata.time,
        **text_common_args,
    )
    ax[1].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
    ax[1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
    ax[1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
    ax[1].text(
        0.5,
        0.4,
        kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
        **text_common_args,
    )
    ax[1].text(
        0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
    )
    ax[1].text(
        0.5,
        0.2,
        "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
        **text_common_args,
    )
    ax[1].text(
        0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
    )
    ax[1].text(
        0.5, 0.0, "Number of particles $N_p$: $%.0f$ " % (Np), **text_common_args
    )
    ax[1].axis("off")

    fig.tight_layout()

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
    # hist_x, _ = np.histogram(data_x, bins=bins_x)
    # hist_y, _ = np.histogram(data_y, bins=bins_y)
    # hist_uncorr = np.outer(hist_x, hist_y).astype("float64")

    # norm
    min_level = 1 / np.sum(hist)
    hist /= np.sum(hist)
    # hist_uncorr /= np.sum(hist_uncorr)

    # print(np.sum(hist),np.sum(hist_uncorr))
    # regularize
    # hist[hist < min_level] = min_level
    # hist_uncorr[hist_uncorr < min_level] = min_level

    # build correlation_funciton
    # correlation_funcion = hist - hist_uncorr

    X, Y = np.meshgrid(xedges, yedges)
    return X, Y, hist  # hist_uncorr, correlation_funcion


def plot_histogram_density_plot(quantity_x, quantity_y, qname_x, qname_y, Nbins=100):
    quantity_x = quantity_x.value
    quantity_y = quantity_y.value

    X, Y, hist = prepare_histogram_density_plot(quantity_x, quantity_y, Nbins)
    nx = 2
    ny = 1
    fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5 * nx, 5 * ny))

    fig.suptitle(f"({qname_x},{qname_y}) space")

    ax[0].set_title(r"$f_{12}(x,y)$")
    ax[0].set_ylabel(qname_y)
    ax[0].set_xlabel(qname_x)
    ax[0].set_box_aspect(1)
    # plt.imshow(hist.T, origin='lower', aspect='auto', cmap='plasma')

    pax = ax[0].pcolormesh(X, Y, hist.T, cmap=prefered_color, norm=LogNorm())
    ax[0].grid(True, linewidth=0.5, alpha=0.3)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")

    x = np.logspace(np.log10(min(quantity_x)), np.log10(max(quantity_x)), 100)
    y = min(quantity_y) / 10 * (x / min(quantity_x)) ** (2 / 3)

    # ax[0].plot(x,y,color='red',alpha=0.5)

    # add panel with infromation about the run
    Np = len(quantity_x)
    text_common_args = dict(
        fontsize=10, ha="center", va="center", transform=ax[1].transAxes
    )

    ax[1].text(
        0.5,
        0.8,
        "Cooling halo with spin at time $t=%.2f$ Gyr" % data.metadata.time,
        **text_common_args,
    )
    ax[1].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
    ax[1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
    ax[1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
    ax[1].text(
        0.5,
        0.4,
        kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
        **text_common_args,
    )
    
    ax[1].text(
        0.5, 0.0, "Number of particles $N_p$: $%.0f$ " % (Np), **text_common_args
    )
    ax[1].axis("off")

    # fig.colorbar(pax)
    fig.tight_layout()
    plt.savefig(sys.argv[2], dpi=200)


if to_plot == "correlation_hist":
    q_x = nH.to(1 / unyt.cm ** 3)
    q_y = T.to(
        unyt.K
    )  # normB.to(1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s)) #T.to(unyt.K)
    q_x = q_x.to_physical()
    q_y = q_y.to_physical()
    plot_histogram_density_plot(
        q_x, q_y, r"$\rho / m_H$ $[1/cm^3]$", "T $[K]$", Nbins=200
    )
    exit()
