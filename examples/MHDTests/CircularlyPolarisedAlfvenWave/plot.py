from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import numpy as np
import sys
import matplotlib.pyplot as plt

snap = int(sys.argv[1])

for dir in sys.argv[2:]:

    file = dir + "CircularlyPolarisedAlfvenWave_%04d.hdf5" % snap
    data = load(file)

    time = data.metadata.time
    eta = data.metadata.parameters["MHD:diffusion_eta"]
    eta = float(eta)

    pos = data.gas.coordinates
    B = data.gas.magnetic_flux_density

    x1 = (pos[:, 0] + 2 * pos[:, 1] + 2 * pos[:, 2]) / 3
    B2 = (B[:, 1] - 2 * B[:, 0]) / np.sqrt(5)

    plt.scatter(x1, B2, s=0.5, label=dir)

k = 2.0 * np.pi
x1exact = np.linspace(0.0, 3.0, 1000)
B2exact = 0.1 * np.sin(2 * np.pi * x1exact) * np.exp(-eta * k * k * time)

plt.plot(x1exact, B2exact, "k-", lw=0.5)
plt.scatter(x1, B2, c="k", s=0.5)
plt.xlabel("x1")
plt.ylabel("B2")

plt.savefig("AlfevWaves", dpi=200)
