import numpy as np
import h5py
import argparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import scipy.stats as stats
import unyt


def bin_vals(x, y, log=True):
    stat = "median" if log else "mean"
    hist, edges, _ = stats.binned_statistic(x, y, statistic=stat, bins=100)
    mids = 0.5 * (edges[1:] + edges[:-1])
    return mids, hist


argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

x = None
rho = None
T = None
v = None
time = None
with h5py.File(args.input, "r") as handle:
    coords = handle["PartType0/Coordinates"][:]
    rho = handle["PartType0/Densities"][:]
    T = handle["PartType0/Temperatures"][:]
    vs = handle["PartType0/Velocities"][:]
    time = handle["Header"].attrs["Time"][0]
    box = handle["Header"].attrs["BoxSize"][:]

    units = dict(handle["Units"].attrs)
    uM = (units["Unit mass in cgs (U_M)"][0] * unyt.g).in_base("galactic")
    uL = (units["Unit length in cgs (U_L)"][0] * unyt.cm).in_base("galactic")
    ut = (units["Unit time in cgs (U_t)"][0] * unyt.s).in_base("galactic")

    coords = (coords * uL).in_base("galactic")
    rho = (rho * uM / uL ** 3).in_base("galactic")
    T = T * unyt.K
    vs = (vs * uL / ut).in_base("galactic")
    time = (time * ut).in_base("galactic")
    box = (box * uL).in_base("galactic")

    coords -= 0.5 * box[None, :]

    x = np.sqrt((coords ** 2).sum(axis=1))
    v = (coords * vs).sum(axis=1) / x

rhohist, edges, _ = stats.binned_statistic(x, rho, statistic="median", bins=100)
mids = 0.5 * (edges[1:] + edges[:-1])
rhohist[np.isnan(rhohist)] = 0.0
imax = np.argmax(rhohist)
xmax = mids[imax]

x.name = "radius"
x.convert_to_units("kpc")
v.name = "radial velocity"
v.convert_to_units("km/s")
rho.name = "density"
rho.convert_to_units("g/cm**3")
T.name = "temperature"

fig, ax = pl.subplots(1, 3, figsize=(8, 4), sharex=True)

with unyt.matplotlib_support:
    ax[0].semilogy(x, rho, ".")
    b, h = bin_vals(x, rho)
    ax[0].semilogy(b, h, "--")

    ax[1].plot(x, v, ".")
    b, h = bin_vals(x, v, log=False)
    ax[1].plot(b, h, "--")

    ax[2].semilogy(x, T, ".")
    b, h = bin_vals(x, T)
    ax[2].semilogy(b, h, "--")

for a in ax:
    a.axvline(x=xmax, color="k", linestyle="--")

ax[0].set_ylim(1.0e-28, 1.0e-24)
ax[1].set_ylim(-10.0, 200.0)
ax[2].set_ylim(10.0, 1.0e8)

ax[1].set_title(f"t = {time:.2e}")

pl.tight_layout()
pl.savefig(args.output, dpi=300)
