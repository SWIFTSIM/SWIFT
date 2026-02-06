################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

# Computes the analytical solution of the Zeldovich pancake and compares with
# the simulation result

# Parameters
T_i = 100.0  # Initial temperature of the gas (in K)
z_c = 1.0  # Redshift of caustic formation (non-linear collapse)
z_i = 100.0  # Initial redshift
gas_gamma = 5.0 / 3.0  # Gas adiabatic index
N_output = 119

# Physical constants needed for internal energy to temperature conversion
kB_in_SI = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

# Read the simulation data
sim = h5py.File("box_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
redshift = sim["/Header"].attrs["Redshift"][0]
a = sim["/Header"].attrs["Scale-factor"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
MHD_scheme = sim["/HydroScheme"].attrs["MHD Scheme"]
git = sim["Code"].attrs["Git Revision"]
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
m_gas = sim["/PartType0/Masses"][0]
N = sim["/Header"].attrs["NumPart_Total"][0]
sim.close()

# Expected comoving quantities
rho_0 = N * m_gas / boxSize ** 3
u_0 = kB_in_SI * T_i / (gas_gamma - 1.0) / mH_in_kg
u_0 *= 1e-6  # conversion to internal units
u_0 *= a ** (-3 * (1 - gas_gamma))
S_0 = (gas_gamma - 1.0) * u_0 * rho_0 ** (-(gas_gamma - 1.0))

# Mean quantities over time
z = np.zeros(N_output)
a = np.zeros(N_output)
S_mean = np.zeros(N_output)
S_std = np.zeros(N_output)
u_mean = np.zeros(N_output)
u_std = np.zeros(N_output)
P_mean = np.zeros(N_output)
P_std = np.zeros(N_output)
rho_mean = np.zeros(N_output)
rho_std = np.zeros(N_output)

vx_mean = np.zeros(N_output)
vy_mean = np.zeros(N_output)
vz_mean = np.zeros(N_output)
vx_std = np.zeros(N_output)
vy_std = np.zeros(N_output)
vz_std = np.zeros(N_output)

bx_mean = np.zeros(N_output)
by_mean = np.zeros(N_output)
bz_mean = np.zeros(N_output)
bx_std = np.zeros(N_output)
by_std = np.zeros(N_output)
bz_std = np.zeros(N_output)

for i in range(N_output):
    sim = h5py.File("box_%04d.hdf5" % i, "r")

    z[i] = sim["/Cosmology"].attrs["Redshift"][0]
    a[i] = sim["/Cosmology"].attrs["Scale-factor"][0]

    S = sim["/PartType0/Entropies"][:]
    S_mean[i] = np.mean(S)
    S_std[i] = np.std(S)

    u = sim["/PartType0/InternalEnergies"][:]
    u_mean[i] = np.mean(u)
    u_std[i] = np.std(u)

    P = sim["/PartType0/Pressures"][:]
    P_mean[i] = np.mean(P)
    P_std[i] = np.std(P)

    rho = sim["/PartType0/Densities"][:]
    rho_mean[i] = np.mean(rho)
    rho_std[i] = np.std(rho)

    v = sim["/PartType0/Velocities"][:, :]
    vx_mean[i] = np.mean(v[:, 0])
    vy_mean[i] = np.mean(v[:, 1])
    vz_mean[i] = np.mean(v[:, 2])
    vx_std[i] = np.std(v[:, 0])
    vy_std[i] = np.std(v[:, 1])
    vz_std[i] = np.std(v[:, 2])

    b = sim["/PartType0/Bfield"][:, :]
    bx_mean[i] = np.mean(b[:, 0])
    by_mean[i] = np.mean(b[:, 1])
    bz_mean[i] = np.mean(b[:, 2])
    bx_std[i] = np.std(b[:, 0])
    by_std[i] = np.std(b[:, 1])
    bz_std[i] = np.std(b[:, 2])

# Move to physical quantities
rho_mean_phys = rho_mean / a ** 3
u_mean_phys = u_mean / a ** (3 * (gas_gamma - 1.0))
S_mean_phys = S_mean
B_mean_phys = (
    sqrt(bx_mean * bx_mean + by_mean * by_mean + bz_mean * bz_mean)
    / sqrt(8 * 3.14)
    / a ** (3.0 / 2.0 * (gas_gamma - 1.0))
)

vx_mean_phys = vx_mean / a
vy_mean_phys = vy_mean / a
vz_mean_phys = vz_mean / a
vx_std_phys = vx_std / a
vy_std_phys = vy_std / a
vz_std_phys = vz_std / a

# Plot the interesting quantities
plt.figure(figsize=(7, 7 / 1.6))

line_color = "C4"
binned_color = "C2"
binned_marker_size = 4

scatter_props = dict(
    marker=".",
    ms=4,
    markeredgecolor="none",
    alpha=0.2,
    zorder=-1,
    rasterized=True,
    linestyle="none",
)

# Density evolution --------------------------------
plt.subplot(231)  # , yscale="log")
plt.semilogx(a, rho_mean / rho_0, **scatter_props)
plt.semilogx([1e-10, 1e10], np.ones(2), "-", color="0.6", lw=1.0)
plt.semilogx([1e-10, 1e10], np.ones(2) * 0.99, "--", color="0.6", lw=1.0)
plt.semilogx([1e-10, 1e10], np.ones(2) * 1.01, "--", color="0.6", lw=1.0)
plt.text(1e-2, 1.0105, "+1%", color="0.6", fontsize=9)
plt.text(1e-2, 0.9895, "-1%", color="0.6", fontsize=9, va="top")
plt.text(1e-2, 1.015, "$\\rho_0=%.3f$" % rho_0)
plt.ylim(0.98, 1.02)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
plt.ylabel("${\\rm Comoving~density}~\\rho / \\rho_0$", labelpad=0.0)

# Thermal energy evolution --------------------------------
plt.subplot(232)  # , yscale="log")
plt.semilogx(a, u_mean / u_0, **scatter_props)
plt.semilogx([1e-10, 1e10], np.ones(2), "-", color="0.6", lw=1.0)
plt.semilogx([1e-10, 1e10], np.ones(2) * 0.99, "--", color="0.6", lw=1.0)
plt.semilogx([1e-10, 1e10], np.ones(2) * 1.01, "--", color="0.6", lw=1.0)
plt.text(1e-2, 1.0105, "+1%", color="0.6", fontsize=9)
plt.text(1e-2, 0.9895, "-1%", color="0.6", fontsize=9, va="top")
plt.text(1e-2, 1.015, "$u_0=%.3e$" % (u_0))
plt.ylim(0.98, 1.02)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
plt.ylabel("${\\rm Comoving~internal~energy}~u / u_0$", labelpad=0.0)

# Entropy evolution --------------------------------
plt.subplot(233)  # , yscale="log")
plt.semilogx(a, S_mean / S_0, **scatter_props)
plt.semilogx([1e-10, 1e10], np.ones(2), "-", color="0.6", lw=1.0)
plt.semilogx([1e-10, 1e10], np.ones(2) * 0.99, "--", color="0.6", lw=1.0)
plt.semilogx([1e-10, 1e10], np.ones(2) * 1.01, "--", color="0.6", lw=1.0)
plt.text(1e-2, 1.0105, "+1%", color="0.6", fontsize=9)
plt.text(1e-2, 0.9895, "-1%", color="0.6", fontsize=9, va="top")
plt.text(1e-2, 1.015, "$A_0=%.3e$" % (S_0))
plt.ylim(0.98, 1.02)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
plt.ylabel("${\\rm Comoving~entropy}~A / A_0$", labelpad=0.0)

# Peculiar velocity evolution ---------------------
plt.subplot(234)
# plt.semilogx(a, vx_mean, **scatter_props)
# plt.semilogx(a, vy_mean, **scatter_props)
# plt.semilogx(a, vz_mean, **scatter_props)
plt.semilogx(a, bx_mean, **scatter_props)
plt.semilogx(a, by_mean, **scatter_props)
plt.semilogx(a, bz_mean, **scatter_props)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
# plt.ylabel("${\\rm Peculiar~velocity~mean}$", labelpad=-5.0)
plt.ylabel("${\\rm Peculiar~B~mean }$", labelpad=-5.0)

# Peculiar velocity evolution ---------------------
plt.subplot(235)
# plt.semilogx(a, vx_std, **scatter_props)
# plt.semilogx(a, vy_std, **scatter_props)
# plt.semilogx(a, vz_std, **scatter_props)
plt.semilogx(a, bx_std, **scatter_props)
plt.semilogx(a, by_std, **scatter_props)
plt.semilogx(a, bz_std, **scatter_props)
# plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
plt.ylabel("${\\rm Peculiar~velocity~std-dev}$", labelpad=0.0)


# Information -------------------------------------
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.plot([-0.45, 0.1], [0.62, 0.62], "k-", lw=1)
plt.text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.45, 0.4, scheme.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.45, 0.3, kernel.decode("utf-8"), fontsize=text_fontsize)
plt.text(
    -0.45,
    0.2,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.xticks([])
plt.yticks([])

plt.tight_layout()

plt.savefig("ConstantBox_comoving.png", dpi=200)


plt.figure(figsize=(7, 7 / 1.6))

# Density evolution --------------------------------
plt.subplot(231)  # , yscale="log")
plt.loglog(a, rho_mean_phys, **scatter_props)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$")
plt.ylabel("${\\rm Physical~density}$")

# Thermal energy evolution --------------------------------
plt.subplot(232)  # , yscale="log")
plt.loglog(a, u_mean_phys, **scatter_props)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$")
plt.ylabel("${\\rm Physical~internal~energy}$")

# Entropy evolution --------------------------------
plt.subplot(233)  # , yscale="log")
plt.semilogx(a, S_mean_phys, **scatter_props)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$")
plt.ylabel("${\\rm Physical~entropy}$")

# Peculiar velocity evolution ---------------------
plt.subplot(234)
plt.semilogx(a, vx_mean_phys, **scatter_props)
plt.semilogx(a, vy_mean_phys, **scatter_props)
plt.semilogx(a, vz_mean_phys, **scatter_props)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
plt.ylabel("${\\rm Peculiar~velocity~mean}$", labelpad=-5.0)

# Peculiar velocity evolution ---------------------
plt.subplot(235)
plt.semilogx(a, vx_std_phys, **scatter_props)
plt.semilogx(a, vy_std_phys, **scatter_props)
plt.semilogx(a, vz_std_phys, **scatter_props)
plt.xlim(8e-3, 1.1)
plt.xlabel("${\\rm Scale-factor}$", labelpad=0.0)
plt.ylabel("${\\rm Peculiar~velocity~std-dev}$", labelpad=0.0)

# Information -------------------------------------
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.plot([-0.45, 0.1], [0.62, 0.62], "k-", lw=1)
plt.text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.45, 0.4, scheme.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.45, 0.3, kernel.decode("utf-8"), fontsize=text_fontsize)
plt.text(
    -0.45,
    0.2,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.xticks([])
plt.yticks([])

plt.tight_layout()

plt.savefig("ConstantBox_physical.png", dpi=200)
