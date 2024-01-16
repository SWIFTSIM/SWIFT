from swiftsimio import load
import numpy as np
import sys
import matplotlib.pyplot as plt

res = [16, 32

x1exact = np.linspace(0.0, 3.0, 1000)
B2exact = 0.1 * np.sin(2 * np.pi * x1exact)
plt.plot(x1exact, B2exact, 'r--', lw=0.5)

for ind, file in enumerate(sys.argv[1:]):

    # Get particle attribuets
    data = load(file)
    B = data.gas.magnetic_flux_densities
    B2 = (B[:, 1] - 2 * B[:, 0]) / np.sqrt(5)

    L1 = 10 * np.sum(abs(B2 - B2exact)) / len(B2)
    plt.scatter(x1, B2, c='k', s=0.5, label=file)

plt.xlabel(r"$x_1$")
plt.ylabel(r"$B_2$")

plt.legend()

plt.savefig("CircularlyPolarisedAlfvenWave.png", dpi=200)
