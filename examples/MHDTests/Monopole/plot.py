from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker

filename = sys.argv[1]
data = load(filename)

divB = data.gas.magnetic_divergences
data.gas.mass_weighted_divB = data.gas.masses * divB

mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)

mass_weighted_divB_map = project_gas(
    data, resolution=1024, project="mass_weighted_divB", parallel=True
)

divB_map = mass_weighted_divB_map / mass_map

from matplotlib.pyplot import imsave

imsave(sys.argv[2], np.rot90(divB_map.value), cmap="gray", vmin=-1.5, vmax=1.5)
