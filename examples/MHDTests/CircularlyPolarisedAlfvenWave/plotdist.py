from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import numpy as np
import sys
import matplotlib.pyplot as plt

filename = sys.argv[1]
data = load(filename)

# Get particle attribuets
pos = data.gas.coordinates
B = data.gas.magnetic_flux_densities

x1 = (pos[:, 0] + 2 * pos[:, 1] + 2 * pos[:, 2]) / 3

B2 = (B[:, 1] - 2 * B[:, 0]) / np.sqrt(5)

x1exact = np.linspace(0.0, 3.0, 1000)
B2exact = 0.1 * np.sin(2 * np.pi * x1exact)

plt.scatter(pos[:, 0], pos[:, 1], c=B2, s=1.0)
# plt.scatter(x1, B2, c='k', s=0.5)
plt.xlabel("x")
plt.ylabel("y")

plt.savefig(sys.argv[2], dpi=200)
