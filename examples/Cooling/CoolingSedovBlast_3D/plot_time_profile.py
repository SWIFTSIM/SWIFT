import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse
import unyt

argparser = argparse.ArgumentParser()
argparser.add_argument("input", nargs="+")
argparser.add_argument("output")
argparser.add_argument("--names", "-n", nargs="+")
args = argparser.parse_args()

fig, ax = pl.subplots(1, 2, sharex=True)

for ifile, file in enumerate(args.input):
    data = np.loadtxt(
        file,
        dtype=[("time", np.float32), ("radius", np.float32), ("velocity", np.float32)],
    )

    t = data["time"] * unyt.Myr
    t.name = "time"
    r = data["radius"] * unyt.kpc
    r.name = "radius"
    v = data["velocity"] * unyt.km / unyt.s
    v.name = "velocity"

    label = None
    if args.names:
        label = args.names[ifile]
    with unyt.matplotlib_support:
        ax[0].plot(t, r, "-", label=label)
        ax[1].plot(t, v, "-", label=label)

if args.names:
    ax[1].legend(loc="best")

pl.tight_layout()
pl.savefig(args.output, dpi=300)
