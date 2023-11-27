import matplotlib.pyplot as plt
import sys

import numpy as np
import swiftsimio
from swiftsimio.visualisation import project_gas
import scipy.interpolate as si
from pathlib import Path


def project_nn(x, y, data, data_xy):
    xv, yv = np.meshgrid(x, y)
    return si.griddata(data_xy, data, (xv, yv), method="nearest")

if __name__ == "__main__":
    try:
        n = sys.argv[1]
    except IndexError:
        n = 2

    cwd = Path(__file__).parent

    fname = cwd / f"output_{n:04}.hdf5"
    # fname = cwd / f"IC.hdf5"

    data = swiftsimio.load(fname)

    # parameters for swiftsimio projections
    projection_kwargs = {"resolution": 1024, "parallel": True}

    mass_map = project_gas(data, project="masses", **projection_kwargs)
    mass_weighted_metal_map = project_gas(data, project="metal_mass_fractions", **projection_kwargs)
    metal_map = mass_weighted_metal_map / mass_map

    res = 1000
    plt.imshow(metal_map.T)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("metal_advection.png", dpi=300)
    plt.show()