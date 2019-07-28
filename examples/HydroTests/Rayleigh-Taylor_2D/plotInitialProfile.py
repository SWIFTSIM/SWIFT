#!/usr/bin/env python
import makeIC
import matplotlib.pyplot as plt
from swiftsimio import load
import numpy as np


N = 1000
filename = "rayleigh_taylor_0000.hdf5"

f = load(filename)

# Get data from snapshot
x, y, _ = f.gas.coordinates.value.T
rho = f.gas.density.value
a = f.gas.entropy.value

# Get analytical solution
y_an = np.linspace(0, makeIC.box_size[1], N)
rho_an = makeIC.density(y_an)
a_an = makeIC.entropy(y_an)

plt.figure()
plt.plot(a, y, ".r", label="Entropy")
plt.plot(a_an, y_an, "--r")

plt.plot(rho, y, ".k", label="Density")
plt.plot(rho_an, y_an, "--k")

plt.ylabel("Position")
plt.xlabel("Density / Entropy")

plt.xlim(0, 2.5)
plt.ylim(0, makeIC.box_size[1])
plt.legend()

plt.show()
