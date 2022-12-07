import numpy as np
import argparse
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

data = np.loadtxt(
    args.input,
    usecols=(1, 13, 14, 16),
    dtype=[
        ("Time", np.float32),
        ("Ekin", np.float32),
        ("Etherm", np.float32),
        ("Erad", np.float32),
    ],
)

pl.plot(data["Time"], data["Ekin"], label="Kinetic energy")
pl.plot(data["Time"], data["Etherm"], label="Thermal energy")
pl.plot(data["Time"], data["Erad"], label="Radiated energy")
pl.plot(
    data["Time"], data["Ekin"] + data["Etherm"] + data["Erad"], label="Total energy"
)

pl.ylabel("Energy")
pl.xlabel("Time")

pl.legend(loc="best")

pl.tight_layout()
pl.savefig(args.output, dpi=300)
