import numpy as np
import argparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

argparser = argparse.ArgumentParser()
argparser.add_argument("input", nargs="+")
argparser.add_argument("output")
argparser.add_argument("--names", "-n", nargs="+")
args = argparser.parse_args()

for ifile, file in enumerate(args.input):

    data = np.loadtxt(
        file,
        usecols=(1, 13, 14, 16),
        dtype=[
            ("Time", np.float32),
            ("Ekin", np.float32),
            ("Etherm", np.float32),
            ("Erad", np.float32),
        ],
    )

    color = f"C{ifile}"

    if args.names:
        pl.plot([], [], "-", color=color, label=args.names[ifile])

    pl.plot(data["Time"], data["Ekin"], "-.", color=color)
    pl.plot(data["Time"], data["Etherm"], "--", color=color)
    pl.plot(data["Time"], data["Erad"], ":", color=color)
    pl.plot(
        data["Time"], data["Ekin"] + data["Etherm"] + data["Erad"], "-", color=color
    )

pl.plot([], [], "k-", label="Total energy")
pl.plot([], [], "k-.", label="Kinetic energy")
pl.plot([], [], "k--", label="Thermal energy")
pl.plot([], [], "k:", label="Radiated energy")
pl.ylabel("Energy")
pl.xlabel("Time")

pl.legend(loc="best", ncol=3)

pl.tight_layout()
pl.savefig(args.output, dpi=300)
