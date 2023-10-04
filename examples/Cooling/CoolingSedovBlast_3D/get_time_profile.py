import numpy as np
import h5py
import argparse
import scipy.stats as stats
import unyt

argparser = argparse.ArgumentParser()
argparser.add_argument("input", nargs="+")
argparser.add_argument("output")
args = argparser.parse_args()

nfile = len(args.input)

time = np.zeros(nfile)
radius = np.zeros(nfile)
velocity = np.zeros(nfile)

for ifile, file in enumerate(args.input):
    with h5py.File(file, "r") as handle:
        coords = handle["PartType0/Coordinates"][:]
        rho = handle["PartType0/Densities"][:]
        vs = handle["PartType0/Velocities"][:]
        t = handle["Header"].attrs["Time"][0]
        box = handle["Header"].attrs["BoxSize"][:]

        units = dict(handle["Units"].attrs)
        uM = (units["Unit mass in cgs (U_M)"][0] * unyt.g).in_base("galactic")
        uL = (units["Unit length in cgs (U_L)"][0] * unyt.cm).in_base("galactic")
        ut = (units["Unit time in cgs (U_t)"][0] * unyt.s).in_base("galactic")

        coords = (coords * uL).in_base("galactic")
        rho = (rho * uM / uL ** 3).in_base("galactic")
        vs = (vs * uL / ut).in_base("galactic")
        t = (t * ut).in_base("galactic")
        box = (box * uL).in_base("galactic")

        coords -= 0.5 * box[None, :]

        x = np.sqrt((coords ** 2).sum(axis=1))
        v = (coords * vs).sum(axis=1) / x

        x.convert_to_units("kpc")
        v.convert_to_units("km/s")
        t.convert_to_units("Myr")

        rhohist, _, _ = stats.binned_statistic(x, rho, statistic="median", bins=100)
        vhist, edges, _ = stats.binned_statistic(x, v, statistic="mean", bins=100)
        mids = 0.5 * (edges[1:] + edges[:-1])

        rhohist[np.isnan(rhohist)] = 0.0
        imax = np.argmax(rhohist)

        time[ifile] = t
        radius[ifile] = mids[imax]
        velocity[ifile] = vhist[imax]

isort = np.argsort(time)

time = time[isort]
radius = radius[isort]
velocity = velocity[isort]

with open(args.output, "w") as handle:
    handle.write("# Radial profile of shock wave\n")
    handle.write("# Column 0: time (Myr)\n")
    handle.write("# Column 1: radius (kpc)\n")
    handle.write("# Column 2: velocity (km/s)\n")
    for t, r, v in zip(time, radius, velocity):
        handle.write(f"{t:.5e}\t{r:.5e}\t{v:.5e}\n")
