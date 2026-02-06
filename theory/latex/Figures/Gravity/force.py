import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

# tex stuff
# rc('font',**{'family':'serif','serif':['Palatino']})
params = {
    "axes.labelsize": 14,
    "font.size": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "text.usetex": True,
    "figure.figsize": (9, 6),
    "figure.subplot.left": 0.08,  # the left side of the subplots of the figure
    "figure.subplot.right": 0.95,  # the right side of the subplots of the figure
    "figure.subplot.bottom": 0.1,  # the bottom of the subplots of the figure
    "figure.subplot.top": 0.93,  # the top of the subplots of the figure
    "figure.subplot.wspace": 0.2,  # the amount of width reserved for blank space between subplots
    "figure.subplot.hspace": 0.2,  # the amount of height reserved for white space between subplots
    "lines.markersize": 4,
    #'axes.formatter.limits' : (, 0),
}
matplotlib.rcParams.update(params)
matplotlib.rc("font", family="serif")
import sys
import math
import os


epsilon = 1.0
h = 2.8 * epsilon
r_s = 20.0
r_c = 4.5 * r_s

r = np.arange(0.1, 2500, 1 / 1000.0)

numElem = r.size
f = np.zeros(numElem)
f_wanted = np.zeros(numElem)
fac = np.zeros(numElem)
f2 = np.zeros(numElem)
fac2 = np.zeros(numElem)
kernel = np.zeros(numElem)

for i in range(numElem):
    if r[i] >= r_c:
        f[i] = 0.0
        f_wanted[i] = 1.0 / r[i] ** 3
    elif r[i] >= h:
        f[i] = 1.0 / r[i] ** 3
        f_wanted[i] = 1.0 / r[i] ** 3
    else:
        u = r[i] / h
        if u < 0.5:
            f[i] = (10.666666666667 + u * u * (32.0 * u - 38.4)) / h ** 3
            f_wanted[i] = (10.666666666667 + u * u * (32.0 * u - 38.4)) / h ** 3
        else:
            f[i] = (
                21.333333333333
                - 48.0 * u
                + 38.4 * u * u
                - 10.666666666667 * u * u * u
                - 0.066666666667 / (u * u * u)
            ) / h ** 3
            f_wanted[i] = (
                21.333333333333
                - 48.0 * u
                + 38.4 * u * u
                - 10.666666666667 * u * u * u
                - 0.066666666667 / (u * u * u)
            ) / h ** 3

    fac[i] = math.erfc(r[i] / (2.0 * r_s)) + (r[i] / (r_s * np.sqrt(np.pi))) * np.exp(
        -r[i] * r[i] / (4 * r_s * r_s)
    )
    fac2[i] = math.erf(r[i] / (2.0 * r_s)) - (r[i] / (r_s * np.sqrt(np.pi))) * np.exp(
        -r[i] * r[i] / (4 * r_s * r_s)
    )
    f[i] = f[i] * fac[i]
    f2[i] = (1.0 / r[i] ** 3) * fac2[i]

for i in range(numElem):
    u = r[i] / h
    if u < 0.5:
        kernel[i] = (10.666666666667 + u * u * (32.0 * u - 38.4)) / h ** 3
    else:
        kernel[i] = (
            21.333333333333
            - 48.0 * u
            + 38.4 * u * u
            - 10.666666666667 * u * u * u
            - 0.066666666667 / (u * u * u)
        ) / h ** 3

plt.figure()
plt.loglog(r / h, 1 / r ** 3, "b-", label="Newton's law")
plt.loglog(r / h, kernel, "r-", label="Softening kernel")
plt.loglog(r / h, fac / r ** 3, "g-", label="Unsoftend tree force")
plt.loglog(r / h, f2, "m-", label="Mesh force")
plt.loglog(r / h, f, "c--", label="Total tree force", linewidth=3)
plt.loglog(r / h, f + f2, "k--", label="Total force", linewidth=3)


plt.plot([1, 1], [1e-10, 1e10], "k--")
plt.text(
    0.85, 5e-9, r"$h_c=2.8\epsilon$", rotation="vertical", fontsize=14, va="bottom"
)

plt.plot([epsilon / h, epsilon / h], [1e-10, 1e10], "k--")
plt.text(
    0.85 * epsilon / h,
    5e-9,
    r"$\epsilon$",
    rotation="vertical",
    fontsize=14,
    va="bottom",
)

plt.plot([r_s / h, r_s / h], [1e-10, 1e10], "k--")
plt.text(0.85 * r_s / h, 5e-4, "$r_s$", rotation="vertical", fontsize=14, va="bottom")

plt.plot([r_c / h, r_c / h], [1e-10, 1e10], "k--")
plt.text(
    0.85 * r_c / h, 5e-4, "$r_c=4.5r_s$", rotation="vertical", fontsize=14, va="bottom"
)


plt.legend(loc="upper right", ncol=2)


plt.grid()

plt.xlim(1e-1, 200)
plt.xlabel("$r/h_c$")
plt.xticks([0.1, 1, 10, 100], ["$0.1$", "$1$", "$10$", "$100$"])

plt.ylim(4e-10, 1e2)
plt.ylabel("Acceleration")


plt.savefig("force.png")


plt.figure()

plt.semilogx(r / r_s, (f + f2) / f_wanted - 1, "r-", label="$error$", linewidth=2)

plt.plot([1, 1], [-1e10, 1e10], "k--")
plt.text(0.85 * 1, 0.011, "$r_s$", rotation="vertical", fontsize=14, va="bottom")

plt.plot([r_c / r_s, r_c / r_s], [-1e10, 1e10], "k--")
plt.text(
    0.85 * r_c / r_s,
    0.011,
    "$r_c=4.5r_s$",
    rotation="vertical",
    fontsize=14,
    va="bottom",
)

plt.grid()

plt.xlim(1e-1, 200)
plt.xlabel("$r/r_s$")
plt.xticks([0.1, 1, 10, 100], ["$0.1$", "$1$", "$10$", "$100$"])

plt.ylim(-0.025, 0.025)
plt.yticks(
    [-0.02, -0.01, 0.0, 0.01, 0.02],
    [r"$-2\%$", r"$-1\%$", r"$0\%$", r"$1\%$", r"$2\%$"],
)

plt.savefig("error.png")


plt.figure()
plt.loglog(r / r_s, fac, "b-", label="$f_{LR}$", linewidth=2)
plt.loglog(r / r_s, 1 - fac, "r-", label="$1-f_{LR}$", linewidth=2)

plt.grid()
plt.xlim(0.05, 30)
plt.ylim(1e-3, 2)

plt.xlabel("$r/r_s$")
plt.xticks([0.1, 1, 10], ["$0.1$", "$1$", "$10$"])

plt.yticks([0.001, 0.01, 0.1, 1], ["$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "$1$"])

plt.legend(loc="lower left")

plt.plot([1, 1], [1e-10, 1e10], "k--")
plt.text(0.85 * 1, 2e-3, "$r_s$", rotation="vertical", fontsize=14, va="bottom")

plt.plot([r_c / r_s, r_c / r_s], [1e-10, 1e10], "k--")
plt.text(
    0.85 * r_c / r_s,
    2e-3,
    "$r_c=4.5r_s$",
    rotation="vertical",
    fontsize=14,
    va="bottom",
)


plt.savefig("correction.png")
