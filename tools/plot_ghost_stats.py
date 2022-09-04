#!/usr/bin/env python
#
# Usage:
#  python plot_ghost_stats.py ghost_stats_file_1 [ghost_stats_file_2...]
#    output_image.png
#
# Description:
#  Plot the ghost statistics contained in the given input files. For each
#  particle type (hydro, stars, black holes), three plots are generated:
#   1. A histogram of the number of active particles as a function of iteration
#      number.
#   2. The distribution of h values for active particles as a function of the
#      iteration number, compared with the final (converged) values.
#   3. The number of particles that are treated as having no neighbours as a
#      function of the iteration number.

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument("input", nargs="+", action="store")
argparser.add_argument("output", action="store")
argparser.add_argument("--label", "-l", action="store", default=None)
args = argparser.parse_args()

data = None
nbin = None
nval = None
for input in args.input:
    # read the header block to figure out how many blocks per particle type
    # we have
    with open(input, "r") as ifile:
        header = ifile.readlines()[:15]
        this_nbin = int(header[3].split()[-1])
        this_nval = int(header[4].split()[-1])
    if nbin is None:
        nbin = this_nbin
    else:
        if nbin != this_nbin:
            print("Number of bins in files does not match!")
            exit(1)
    if nval is None:
        nval = this_nval
    else:
        if nval != this_nval:
            print("Number of values per bin in files does not match!")
            exit(1)

    tdata = np.loadtxt(input, ndmin=2)
    if len(tdata) > 0:
        if data is None:
            data = np.array(tdata)
        else:
            data = np.append(data, tdata, axis=0)

# omit meaningless first column
data = data[:, 1:]

hdata = data[:, : nval * nbin]
sdata = data[:, nval * nbin : 2 * nval * nbin]
bdata = data[:, 2 * nval * nbin :]


def get_total_stats(pdata):
    counts = pdata[:, ::nval]
    ncell = (counts[:, 0] > 0).sum()
    if ncell == 0:
        return None, None, None, None, None, None, None, 0
    nongb_counts = np.sum(pdata[:, 1::nval], axis=0)
    hmin = np.min(pdata[:, 2::nval], axis=0, initial=np.inf, where=counts > 0)
    hmax = np.max(pdata[:, 3::nval], axis=0, initial=0.0, where=counts > 0)
    hsum = np.sum(pdata[:, 4::nval], axis=0, where=counts > 0)
    hsum2 = np.sum(pdata[:, 5::nval], axis=0, where=counts > 0)
    counts = counts.sum(axis=0)
    filt = counts > 0
    hmean = hsum
    hmean[filt] /= counts[filt]
    hmean2 = hsum2
    hmean2[filt] /= counts[filt]
    hstd = np.sqrt(hmean2 - hmean ** 2)
    bins = np.arange(hmean.shape[0])
    return (
        bins[filt],
        counts[filt],
        nongb_counts[filt],
        hmin[filt],
        hmax[filt],
        hmean[filt],
        hstd[filt],
        ncell,
    )


def no_values(ax):
    ax.text(0.5, 0.5, "No values", horizontalalignment="center", transform=ax.transAxes)
    return


hbins, hcounts, hno_ngb, hhmin, hhmax, hhmean, hhstd, hncell = get_total_stats(
    hdata[:, : nval * (nbin - 1)]
)
htbins, htcounts, htno_ngb, hthmin, hthmax, hthmean, hthstd, htncell = get_total_stats(
    hdata[:, nval * (nbin - 1) :]
)
sbins, scounts, sno_ngb, shmin, shmax, shmean, shstd, sncell = get_total_stats(
    sdata[:, : nval * (nbin - 1)]
)
stbins, stcounts, stno_ngb, sthmin, sthmax, sthmean, sthstd, stncell = get_total_stats(
    sdata[:, nval * (nbin - 1) :]
)
bbins, bcounts, bno_ngb, bhmin, bhmax, bhmean, bhstd, bncell = get_total_stats(
    bdata[:, : nval * (nbin - 1)]
)
btbins, btcounts, btno_ngb, bthmin, bthmax, bthmean, bthstd, btncell = get_total_stats(
    bdata[:, nval * (nbin - 1) :]
)

fig, ax = pl.subplots(3, 3, sharex=True, figsize=(10, 8))

if hncell > 0:
    ax[0][0].axhline(y=htcounts[0], linestyle="--", color="k")
    ax[1][0].axhline(y=hthmin, linestyle=":", color="k")
    ax[1][0].axhline(y=hthmean - hthstd, linestyle="--", color="k")
    ax[1][0].axhline(y=hthmean, linestyle="-", color="k")
    ax[1][0].axhline(y=hthmean + hthstd, linestyle="--", color="k")
    ax[1][0].axhline(y=hthmax, linestyle=":", color="k")
    ax[2][0].axhline(y=htno_ngb[0], linestyle="--", color="k")

    ax[0][0].bar(hbins, hcounts, log=True)

    ax[1][0].plot(hbins, hhmin, label="min")
    ax[1][0].plot(hbins, hhmax, label="max")
    ax[1][0].fill_between(hbins, hhmean - hhstd, hhmean + hhstd, color="C2", alpha=0.5)
    ax[1][0].plot(hbins, hhmean, color="C2", label="mean")
    ax[1][0].set_yscale("log")

    if hno_ngb.sum() > 0:
        ax[2][0].bar(hbins, hno_ngb, log=True)
    else:
        no_values(ax[2][0])
else:
    no_values(ax[0][0])
    no_values(ax[1][0])
    no_values(ax[2][0])

if sncell > 0:
    ax[0][1].axhline(y=stcounts[0], linestyle="--", color="k")
    ax[1][1].axhline(y=sthmin, linestyle=":", color="k")
    ax[1][1].axhline(y=sthmean - sthstd, linestyle="--", color="k")
    ax[1][1].axhline(y=sthmean, linestyle="-", color="k")
    ax[1][1].axhline(y=sthmean + sthstd, linestyle="--", color="k")
    ax[1][1].axhline(y=sthmax, linestyle=":", color="k")
    ax[2][1].axhline(y=stno_ngb[0], linestyle="--", color="k")

    ax[0][1].bar(sbins, scounts, log=True)

    ax[1][1].plot(sbins, shmin)
    ax[1][1].plot(sbins, shmax)
    ax[1][1].fill_between(sbins, shmean - shstd, shmean + shstd, color="C2", alpha=0.5)
    ax[1][1].plot(sbins, shmean, color="C2")
    ax[1][1].set_yscale("log")

    if sno_ngb.sum() > 0:
        ax[2][1].bar(sbins, sno_ngb, log=True)
    else:
        no_values(ax[2][1])
else:
    no_values(ax[0][1])
    no_values(ax[1][1])
    no_values(ax[2][1])

if bncell > 0:
    ax[0][2].axhline(y=btcounts[0], linestyle="--", color="k")
    ax[1][2].axhline(y=bthmin, linestyle=":", color="k")
    ax[1][2].axhline(y=bthmean - bthstd, linestyle="--", color="k")
    ax[1][2].axhline(y=bthmean, linestyle="-", color="k")
    ax[1][2].axhline(y=bthmean + bthstd, linestyle="--", color="k")
    ax[1][2].axhline(y=bthmax, linestyle=":", color="k")
    ax[2][2].axhline(y=btno_ngb[0], linestyle="--", color="k")

    ax[0][2].bar(bbins, bcounts, log=True)

    ax[1][2].plot(bbins, bhmin)
    ax[1][2].plot(bbins, bhmax)
    ax[1][2].fill_between(bbins, bhmean - bhstd, bhmean + bhstd, color="C2", alpha=0.5)
    ax[1][2].plot(bbins, bhmean, color="C2")
    ax[1][2].set_yscale("log")

    if bno_ngb.sum() > 0:
        ax[2][2].bar(bbins, bno_ngb, log=True)
    else:
        no_values(ax[2][2])
else:
    no_values(ax[0][2])
    no_values(ax[1][2])
    no_values(ax[2][2])

ax[1][0].set_xlabel("iteration")
ax[1][1].set_xlabel("iteration")
ax[1][2].set_xlabel("iteration")

ax[0][0].set_ylabel("count")
ax[1][0].set_ylabel("h")
ax[2][0].set_ylabel("no ngb count")

ax[1][0].legend(loc="best")
ax[1][1].plot([], [], "k-", label="converged value")
ax[1][1].legend(loc="best")

ax[0][0].set_xlim(-0.5, nbin - 1.5)
ax[0][0].set_xticks(np.arange(nbin))
ax[0][0].set_title("hydro")
ax[0][1].set_title("stars")
ax[0][2].set_title("black holes")
ax[1][0].set_title("{0} cells".format(hncell))
ax[1][1].set_title("{0} cells".format(sncell))
ax[1][2].set_title("{0} cells".format(bncell))

if not args.label is None:
    pl.suptitle(args.label)
else:
    pl.suptitle(",".join(args.input))

pl.tight_layout()
pl.savefig(args.output, dpi=300)
