#!/usr/bin/env python3

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

matplotlib.style.use("../../../tools/stylesheets/mnras.mplstyle")


# read in file

evals = np.loadtxt("eigenvals.txt")
n_points = int(np.sqrt(evals.shape[0]) + 0.5)
evals = evals.reshape((n_points, n_points, 4))

# get array for thetas and f's
f = np.linspace(0, 1, num=n_points, endpoint=True)
theta = np.linspace(0, 1, num=n_points, endpoint=True)

f1ind = n_points - 1
theta0ind = 0
thetapihalfind = np.argmin(np.abs(theta - 0.5))


def plot_eigenvals(ax, abscissa, eigenvals):

    ax.plot(abscissa, eigenvals[:, 0], "-", label="$\lambda_1$")
    ax.plot(abscissa, eigenvals[:, 1], "-.", label="$\lambda_2$")
    ax.plot(abscissa, eigenvals[:, 2], "--", label="$\lambda_3$")
    ax.plot(abscissa, eigenvals[:, 3], ":", label="$\lambda_4$")
    return


fig = plt.figure(figsize=(7, 4))
ax1 = fig.add_subplot(1, 3, 1)
plot_eigenvals(ax1, f, evals[:, theta0ind, :])
ax1.set_xlabel("$f$", usetex=True)
ax1.set_title("$\\theta = 0$", usetex=True)
ax1.set_ylabel("Eigenvalues", usetex=True)

ax2 = fig.add_subplot(1, 3, 2)
plot_eigenvals(ax2, f, evals[:, thetapihalfind, :])
ax2.set_xlabel("$f$", usetex=True)
ax2.set_title("$\\theta = \pi/2$", usetex=True)

ax3 = fig.add_subplot(1, 3, 3)
plot_eigenvals(ax3, theta, evals[f1ind, :, :])
ax3.set_xlabel("$\\theta / \pi$", usetex=True)
ax3.set_title("$f = 1$", usetex=True)

for ax in fig.axes:
    ax.set_ylim(-1.01, 1.01)
    ax.set_xlim(-0.01, 1.01)

#  ax3.legend(loc="center left", bbox_to_anchor=(1., 0.5))
ax3.legend()

#  plt.show()
plt.tight_layout()
plt.savefig("HLL_eigenvalues.pdf")
