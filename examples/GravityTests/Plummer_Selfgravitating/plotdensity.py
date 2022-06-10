import numpy as np
from matplotlib import pyplot as plt
import argparse
from scipy.optimize import curve_fit
import h5py

# Parse user input
parser = argparse.ArgumentParser(
    description="Plot multiple density profiles against theoretical prediction"
)
parser.add_argument("files", nargs="+", help="snapshot files to be imaged")
parser.add_argument("--notex", action="store_true", help="Flag to not use LaTeX markup")
parser.add_argument(
    "-a", type=float, default=0.1, help="Softening length of theoretical model"
)
parser.add_argument(
    "-M", type=float, default=1.0e-5, help="Total Mass of theoretical model"
)
parser.add_argument("-Rmin", type=float, default=0.01, help="Min Radius")
parser.add_argument("-Rmax", type=float, default=2.1, help="Max Radius")
parser.add_argument(
    "-nbsamples", type=int, default=64, help="Number of radii to sample (bins)"
)
parser.add_argument(
    "-shift", type=float, default=2.0, help="Shift applied to particles in params.yml"
)

args = parser.parse_args()
fnames = args.files

# Limit number of snapshots to plot
if len(fnames) > 20:
    raise ValueError(
        "Too many ({:d}) files provided (cannot plot more than 20).".format(len(fnames))
    )

# Set parameters
tex = not args.notex
if tex:
    plt.rcParams.update({"text.usetex": tex, "font.size": 14, "font.family": "serif"})
else:
    plt.rcParams.update({"font.size": 12})
figsize = 7

# Model Parameters (Mestel surface density)
G = 4.299581e04
rsp = np.logspace(np.log10(args.Rmin), np.log10(args.Rmax), args.nbsamples)


def plummer_analytical(r):
    return (
        3.0
        * args.M
        / (4.0 * np.pi * args.a ** 3)
        * (1.0 + r ** 2 / args.a ** 2) ** (-2.5)
    )


# Plot densities
fig, ax = plt.subplots(figsize=(1.2 * figsize, figsize))
for fname in fnames:
    print(fname)
    f = h5py.File(fname, "r")
    pos = np.array(f["DMParticles"]["Coordinates"]) - args.shift
    time = (
        f["Header"].attrs["Time"][0]
        * f["Units"].attrs["Unit time in cgs (U_t)"][0]
        / 31557600.0e9
    )
    mass = np.array(f["DMParticles"]["Masses"])
    x = pos[:, 0]
    y = pos[:, 1]
    r = np.sqrt(np.sum(pos ** 2, 1))

    # Methods to compute density profile
    def mass_ins(R):
        return ((r < R) * mass).sum()

    mass_ins_vect = np.vectorize(mass_ins)

    def density(R):
        return np.diff(mass_ins_vect(R)) / np.diff(R) / (4.0 * np.pi * R[1:] ** 2)

    dens = density(rsp)
    rs = rsp[1:]

    # remove empty bins
    c = dens > 0
    dens = np.compress(c, dens)
    rs = np.compress(c, rs)

    # Plot
    # ax.loglog(rsp[1:], density(rsp), "o", ms=1.7, label=r"$t=$ {:.3f} Gyr".format(time))
    ax.plot(rs, dens, label=r"$t=$ {:.3f} Gyr".format(time))

ax.plot(rsp, plummer_analytical(rsp), c="black", label="Analytical Prediction")
ax.set_xlabel("r [kpc]")
ax.legend()
ax.loglog()
ax.set_title(
    r"Plummer Density Profile: $a = {:.1e}$ kpc, $M = {:.1e}$ M$_{{\odot}}$".format(
        args.a, args.M * 1e10
    )
)
plt.tight_layout()
ax.set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]")
plt.savefig("density_profile.png")
