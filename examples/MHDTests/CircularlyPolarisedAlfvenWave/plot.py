from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

# Test parameters
wl = 1.0
k = 2.0 * np.pi / wl

v0 = 0.1
B0 = 0.1

mu0=1

sina = 2.0 / 3.0
cosa = np.sqrt(1 - sina ** 2)

sinb = 2.0 / np.sqrt(5)
cosb = np.sqrt(1 - sinb ** 2)

# Rotation matrix to rotate frame S to the projected direction frame (direction of wave propagation) Sp ...
Rotation_Sp_to_S = np.array(
    [
        [cosa * cosb, -sinb, -sina * cosb],
        [cosa * sinb, cosb, -sina * sinb],
        [sina, 0, cosa],
    ]
)

# ... and its inverse to go from Sp to S
Rotation_S_to_Sp = Rotation_Sp_to_S.T

# File to read
filename = sys.argv[1]

# Read in snapshot
with h5py.File(filename, "r") as handle:
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"]

    rho = handle["PartType0/Densities"][:]

    pos = handle["PartType0/Coordinates"][:]
    v = handle["PartType0/Velocities"][:]

    h = handle["PartType0/SmoothingLengths"][:]

    B = handle["PartType0/MagneticFluxDensities"][:]
    divB = handle["PartType0/MagneticDivergences"][:]
    J = handle["PartType0/MagneticFluxCurl"][:]
    Fl = handle["PartType0/LorentzForces"][:]
# Plot things
fig, ax = plt.subplots(2, 3, figsize=(11, 6))

ax00b = ax[0, 0].twinx()
ax10b = ax[1, 0].twinx()

normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
err = np.log10(h * abs(divB) / normB)

pos = pos.T
v = v.T
B = B.T
J = J.T
Fl = Fl.T

posp = np.dot(Rotation_S_to_Sp, pos)
vp = np.dot(Rotation_S_to_Sp, v)
Bp = np.dot(Rotation_S_to_Sp, B)
Jp = np.dot(Rotation_S_to_Sp, J)
Flp = np.dot(Rotation_S_to_Sp, Fl)


for ind, axi in enumerate(ax[0, :]):
    axi.scatter(posp[0, :], Bp[ind, :], s=0.2, label=filename)

for ind, axi in enumerate(ax[1, :]):
    axi.scatter(posp[0, :], vp[ind, :], s=0.2, label=filename)

ax00b.scatter(posp[0, :], rho, s=0.2, c="g")
ax10b.scatter(posp[0, :], err, s=0.2, c="r")

pts = 1000

x1exact = np.linspace(0.0, max(posp[0, :]), pts)
vexact = [np.zeros(pts), v0 * np.sin(k * x1exact), v0 * np.cos(k * x1exact)]
Bexact = [np.ones(pts), B0 * np.sin(k * x1exact), B0 * np.cos(k * x1exact)]

Jexact = [np.zeros(pts), B0 * k * np.sin(k * x1exact), B0 * k * np.cos(k * x1exact)] 
Flexact = [np.zeros(pts), B0 * k * np.cos(k * x1exact), - B0 * k * np.sin(k * x1exact)]

Mexact = [
[np.ones(pts)*(-1.0+B0**2)/mu0,                                       -B0*np.sin(k * x1exact)/mu0,                                       -B0 * np.cos(k * x1exact)/mu0],
[  -B0*np.sin(k * x1exact)/mu0,(np.ones(pts)+B0**2*(np.cos(k * x1exact)-np.sin(k * x1exact)))/mu0,                            -B0**2 * np.sin(2 * k * x1exact)/(2*mu0)],
[-B0 * np.cos(k * x1exact)/mu0,                          -B0**2 * np.sin(2 * k * x1exact)/(2*mu0), (np.ones(pts)+B0**2*(np.sin(k * x1exact)-np.cos(k * x1exact)))/mu0]
]

for ind, axi in enumerate(ax[0, :]):
    ax[0, ind].plot(x1exact, Bexact[ind], "k-", lw=0.5, label="Exact Solution")
    ax[1, ind].plot(x1exact, vexact[ind], "k-", lw=0.5, label="Exact Solution")

labelsv = [r"$v_1$", r"$v_2$", r"$v_3$"]
labelsB = [r"$B_1$", r"$B_2$", r"$B_3$"]

for ind, axi in enumerate(ax[0, :]):
    ax[1, ind].set_xlabel(r"$x_1$")
    ax[0, ind].set_ylabel(labelsB[ind])
    ax[1, ind].set_ylabel(labelsv[ind])
    ax[0, ind].set_xticks([])

ax00b.set_ylabel(r"$\rho$")
ax10b.set_ylabel(r"$err_{\nabla \cdot B}$")

plt.tight_layout()

plt.savefig("AlfvenWaves.png", dpi=100)


# Plot M
#fig, axs = plt.subplots(3, 3, figsize=(3*8, 3*4), sharex=True)
#fig.subplots_adjust(hspace=0.1,wspace=0.5)




#axs[0, 0].plot(x1exact, Mexact[0][0], "k-", lw=0.5, label="Exact Solution")
#axs[0, 1].plot(x1exact, Mexact[0][1], "k-", lw=0.5, label="Exact Solution")
#axs[0, 2].plot(x1exact, Mexact[0][2], "k-", lw=0.5, label="Exact Solution")
#axs[1, 0].plot(x1exact, Mexact[1][0], "k-", lw=0.5, label="Exact Solution")
#axs[1, 1].plot(x1exact, Mexact[1][1], "k-", lw=0.5, label="Exact Solution")
#axs[1, 2].plot(x1exact, Mexact[1][2], "k-", lw=0.5, label="Exact Solution")
#axs[2, 0].plot(x1exact, Mexact[2][0], "k-", lw=0.5, label="Exact Solution")
#axs[2, 1].plot(x1exact, Mexact[2][1], "k-", lw=0.5, label="Exact Solution")
#axs[2, 2].plot(x1exact, Mexact[2][2], "k-", lw=0.5, label="Exact Solution")


#axs[0, 0].set_ylabel(r'$\rm M_{1,1}$')
#axs[0, 1].set_ylabel(r'$\rm M_{1,2}$')
#axs[0, 2].set_ylabel(r'$\rm M_{1,3}$')
#axs[1, 0].set_ylabel(r'$\rm M_{2,1}$')
#axs[1, 1].set_ylabel(r'$\rm M_{2,2}$')
#axs[1, 2].set_ylabel(r'$\rm M_{2,3}$')
#axs[2, 0].set_ylabel(r'$\rm M_{3,1}$')
#axs[2, 1].set_ylabel(r'$\rm M_{3,2}$')
#axs[2, 2].set_ylabel(r'$\rm M_{3,3}$')

#axs[2, 0].set_xlabel(r"$\rm x_1$")
#axs[2, 1].set_xlabel(r"$\rm x_1$")
#axs[2, 2].set_xlabel(r"$\rm x_1$")


#fig.tight_layout()
#plt.savefig("AlfvenWaves_M.png", dpi=100)



fig, axs = plt.subplots(3, 3, figsize=(3*8, 3*4), sharex=True)

axs[0,0].scatter(posp[0, :], Bp[0, :], s=0.2, label=filename)
axs[0,1].scatter(posp[0, :], Bp[1, :], s=0.2, label=filename)
axs[0,2].scatter(posp[0, :], Bp[2, :], s=0.2, label=filename)

axs[1,0].scatter(posp[0, :], Jp[0, :], s=0.2, label=filename)
axs[1,1].scatter(posp[0, :], Jp[1, :], s=0.2, label=filename)
axs[1,2].scatter(posp[0, :], Jp[2, :], s=0.2, label=filename)

axs[2,0].scatter(posp[0, :], Flp[0, :], s=0.2, label=filename)
axs[2,1].scatter(posp[0, :], Flp[1, :], s=0.2, label=filename)
axs[2,2].scatter(posp[0, :], Flp[2, :], s=0.2, label=filename)

axs[0, 0].plot(x1exact, Bexact[0], "k-", lw=0.5, label="Exact Solution")
axs[0, 1].plot(x1exact, Bexact[1], "k-", lw=0.5, label="Exact Solution")
axs[0, 2].plot(x1exact, Bexact[2], "k-", lw=0.5, label="Exact Solution")
axs[1, 0].plot(x1exact, Jexact[0], "k-", lw=0.5, label="Exact Solution")
axs[1, 1].plot(x1exact, Jexact[1], "k-", lw=0.5, label="Exact Solution")
axs[1, 2].plot(x1exact, Jexact[2], "k-", lw=0.5, label="Exact Solution")
axs[2, 0].plot(x1exact, Flexact[0], "k-", lw=0.5, label="Exact Solution")
axs[2, 1].plot(x1exact, Flexact[1], "k-", lw=0.5, label="Exact Solution")
axs[2, 2].plot(x1exact, Flexact[2], "k-", lw=0.5, label="Exact Solution")


axs[0, 0].set_ylabel(r'$\rm B_{1}$')
axs[0, 1].set_ylabel(r'$\rm B_{2}$')
axs[0, 2].set_ylabel(r'$\rm B_{3}$')
axs[1, 0].set_ylabel(r'$\rm J_{1}$')
axs[1, 1].set_ylabel(r'$\rm J_{2}$')
axs[1, 2].set_ylabel(r'$\rm J_{3}$')
axs[2, 0].set_ylabel(r'$\rm F_{L,1}$')
axs[2, 1].set_ylabel(r'$\rm F_{L,2}$')
axs[2, 2].set_ylabel(r'$\rm F_{L,3}$')

axs[2, 0].set_xlabel(r"$\rm x_1$")
axs[2, 1].set_xlabel(r"$\rm x_1$")
axs[2, 2].set_xlabel(r"$\rm x_1$")

fig.tight_layout()

plt.savefig("AlfvenWaves_B_J_F.png", dpi=100)
