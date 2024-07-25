### This file is part of SWIFT.
# Computes the analytical solution of the Zeldovich pancake and compares with
# the simulation result, for a 3D collapse, with magnetic field

import sys
snap = int(sys.argv[1])
# You can input a different directory where a snapshot is saved; if you don't, use the current one
try:
    path = str(sys.argv[2])
except:
    path = ''
filepath = f"{path}zeldovichPancake_%04d.hdf5" % snap

# Parameters
B_0 = 1.0e-9 # Initial magnetic field strength (in T)
## Redshift of caustic formation (non-linear collapse), now 3 values possible
z_c_x = 5.0
z_c_y = 1.0
z_c_z = 0.5
z_i = 100.0  # Initial redshift
a_c_x = 1/(1+z_c_x)
a_c_y = 1/(1+z_c_y)
a_c_z = 1/(1+z_c_z)
a_i = 1/(1+z_i)
T_i = 100.0  # Initial temperature of the gas (in K)
numPart_1D = 32  # Number of particles along each dimension

# Physical constants needed for internal energy to temperature conversion
k_in_J_K = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os.path

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

# Read the simulation data
sim = h5py.File(filepath, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
redshift = sim["/Header"].attrs["Redshift"][0]
a = sim["/Header"].attrs["Scale-factor"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
alpha = sim["/HydroScheme"].attrs["Alpha viscosity"]
git = sim["Code"].attrs["Git Revision"]
mag_eta = sim["/HydroScheme"].attrs["Resistive Eta"][0]

# Cosmological parameters
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
gas_gamma = sim["/HydroScheme"].attrs["Adiabatic index"][0]

x = sim["/PartType0/Coordinates"][:, 0]
y = sim["/PartType0/Coordinates"][:, 1]
z = sim["/PartType0/Coordinates"][:, 2]
v_x = sim["/PartType0/Velocities"][:, 0]
v_y = sim["/PartType0/Velocities"][:, 1]
v_z = sim["/PartType0/Velocities"][:, 2]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]
m = sim["/PartType0/Masses"][:]
x -= 0.5 * boxSize
y -= 0.5 * boxSize
z -= 0.5 * boxSize

# Also read in the magnetic field
Bx = sim["/PartType0/MagneticFluxDensities"][:, 0]
By = sim["/PartType0/MagneticFluxDensities"][:, 1]
Bz = sim["/PartType0/MagneticFluxDensities"][:, 2]

# Derived parameters
rho_0 = m.sum() / (boxSize ** 3)  # critical density of the box
lambda_i = boxSize  # wavelength of the perturbation

unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
unit_current_in_cgs = sim["/Units"].attrs["Unit current in cgs (U_I)"]

unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs
unit_current_in_si = unit_current_in_cgs
unit_B_in_si = unit_mass_in_si * unit_time_in_si ** (-2) * unit_current_in_si ** (-1)
unit_rho_in_si = unit_mass_in_si * unit_length_in_si ** (-3)

# Solutions
q = np.linspace(-0.5 * lambda_i, 0.5 * lambda_i, 256)  # use this for q_x, q_y, q_z
## Set frequencies for all three directions; the same for now
k_x = 2.0 * np.pi / lambda_i
k_y = 2.0 * np.pi / lambda_i
k_z = 2.0 * np.pi / lambda_i

zfac_x = (1.0 + z_c_x) / (1.0 + redshift)
zfac_y = (1.0 + z_c_y) / (1.0 + redshift)
zfac_z = (1.0 + z_c_z) / (1.0 + redshift)

x_s = q - zfac_x * np.sin(k_x * q) / k_x
y_s = q - zfac_y * np.sin(k_y * q) / k_y
z_s = q - zfac_z * np.sin(k_z * q) / k_z

# For the solution of the density, there is some sublety in how to plot, as it depends on qx, qy and qz now.
# Each particle at a certain x will have different values of density because it has different y and z positions 
# where the density changes now.
# So we will make an upper and lower limit of the plot, that should have all points inbetween.
delta_x = boxSize / numPart_1D  # = delta_y = delta_z
q_max = 0.5*delta_x
q_min = -0.5*lambda_i + 0.5*delta_x
rho_s_x_upper = rho_0 / (1.0 - zfac_x * np.cos(k_x * q)) / (1.0 - zfac_y * np.cos(k_y * q_max)) / (1.0 - zfac_z * np.cos(k_z * q_max))
rho_s_x_lower = rho_0 / (1.0 - zfac_x * np.cos(k_x * q)) / (1.0 - zfac_y * np.cos(k_y * q_min)) / (1.0 - zfac_z * np.cos(k_z * q_min))
rho_s_y_upper = rho_0 / (1.0 - zfac_x * np.cos(k_x * q_max)) / (1.0 - zfac_y * np.cos(k_y * q)) / (1.0 - zfac_z * np.cos(k_z * q_max))
rho_s_y_lower = rho_0 / (1.0 - zfac_x * np.cos(k_x * q_min)) / (1.0 - zfac_y * np.cos(k_y * q)) / (1.0 - zfac_z * np.cos(k_z * q_min))
rho_s_z_upper = rho_0 / (1.0 - zfac_x * np.cos(k_x * q_max)) / (1.0 - zfac_y * np.cos(k_y * q_max)) / (1.0 - zfac_z * np.cos(k_z * q))
rho_s_z_lower = rho_0 / (1.0 - zfac_x * np.cos(k_x * q_min)) / (1.0 - zfac_y * np.cos(k_y * q_min)) / (1.0 - zfac_z * np.cos(k_z * q))

v_x_s = -H_0 * (1.0 + z_c_x) / np.sqrt(1.0 + redshift) * np.sin(k_x * q) / k_x
v_y_s = -H_0 * (1.0 + z_c_y) / np.sqrt(1.0 + redshift) * np.sin(k_y * q) / k_y
v_z_s = -H_0 * (1.0 + z_c_z) / np.sqrt(1.0 + redshift) * np.sin(k_z * q) / k_z

# We will plot T and P only as a function of x
T_s_x_upper = T_i * (((1.0 + redshift) / (1.0 + z_i)) ** 3.0 * rho_s_x_upper / rho_0) ** (2.0 / 3.0)
T_s_x_lower = T_i * (((1.0 + redshift) / (1.0 + z_i)) ** 3.0 * rho_s_x_lower / rho_0) ** (2.0 / 3.0)
P_s_x_upper = T_s_x_upper * rho_s_x_upper*unit_rho_in_si * k_in_J_K / mH_in_kg / (a/a_i)**3  ## rho_s is not SI units! So multiply it
P_s_x_lower = T_s_x_lower * rho_s_x_lower*unit_rho_in_si * k_in_J_K / mH_in_kg / (a/a_i)**3

# For magnetic field
y_term = (1.0 - (1.0 + z_c_y) / (1.0 + z_i) * np.cos(k_y * q)) / (1.0 - (1.0 + z_c_y) / (1.0 + redshift) * np.cos(k_y * q))
z_term = (1.0 - (1.0 + z_c_z) / (1.0 + z_i) * np.cos(k_z * q)) / (1.0 - (1.0 + z_c_z) / (1.0 + redshift) * np.cos(k_z * q))
y_term_max = (1.0 - (1.0 + z_c_y) / (1.0 + z_i) * np.cos(k_y * q_max)) / (1.0 - (1.0 + z_c_y) / (1.0 + redshift) * np.cos(k_y * q_max))
z_term_max = (1.0 - (1.0 + z_c_z) / (1.0 + z_i) * np.cos(k_z * q_max)) / (1.0 - (1.0 + z_c_z) / (1.0 + redshift) * np.cos(k_z * q_max))
y_term_min = (1.0 - (1.0 + z_c_y) / (1.0 + z_i) * np.cos(k_y * q_min)) / (1.0 - (1.0 + z_c_y) / (1.0 + redshift) * np.cos(k_y * q_min))
z_term_min = (1.0 - (1.0 + z_c_z) / (1.0 + z_i) * np.cos(k_z * q_min)) / (1.0 - (1.0 + z_c_z) / (1.0 + redshift) * np.cos(k_z * q_min))

Bx_s_x_upper = B_0 * y_term_max * z_term_max * np.ones_like(q)
Bx_s_x_lower = B_0 * y_term_min * z_term_min * np.ones_like(q)
Bx_s_y_upper = B_0 * y_term * z_term_max
Bx_s_y_lower = B_0 * y_term * z_term_min
Bx_s_z_upper = B_0 * y_term_max * z_term
Bx_s_z_lower = B_0 * y_term_min * z_term


# Plot the interesting quantities
plt.figure(figsize=(9, 9 / 1.6))

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


# Density profiles --------------------------------
plt.subplot(4,3,1)  # , yscale="log")
plt.plot(x, rho / rho_0, **scatter_props)
plt.plot(x_s, rho_s_x_upper / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(x_s, rho_s_x_lower / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Comoving~density}}~\\rho / \\rho_0$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

plt.subplot(4,3,2)  # , yscale="log")
plt.plot(y, rho / rho_0, **scatter_props)
plt.plot(y_s, rho_s_y_upper / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(y_s, rho_s_y_lower / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Comoving~density}}~\\rho / \\rho_0$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~y~{\\rm{(Mpc)}}$", labelpad=0)

plt.subplot(4,3,3)  # , yscale="log")
plt.plot(z, rho / rho_0, **scatter_props)
plt.plot(z_s, rho_s_z_upper / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(z_s, rho_s_z_lower / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Comoving~density}}~\\rho / \\rho_0$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~z~{\\rm{(Mpc)}}$", labelpad=0)

# Velocity profiles --------------------------------
plt.subplot(4,3,4)
plt.plot(x, v_x, **scatter_props)
plt.plot(x_s, v_x_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)
plt.ylabel("${\\rm{Peculiar~velocity}}~v_x~{\\rm{(km/s)}}$", labelpad=0)

plt.subplot(4,3,5)
plt.plot(y, v_y, **scatter_props)
plt.plot(y_s, v_y_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Comoving~position}}~y~{\\rm{(Mpc)}}$", labelpad=0)
plt.ylabel("${\\rm{Peculiar~velocity}}~v_y~{\\rm{(km/s)}}$", labelpad=0)

plt.subplot(4,3,6)
plt.plot(z, v_z, **scatter_props)
plt.plot(z_s, v_z_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Comoving~position}}~z~{\\rm{(Mpc)}}$", labelpad=0)
plt.ylabel("${\\rm{Peculiar~velocity}}~v_z~{\\rm{(km/s)}}$", labelpad=0)

# Temperature profile -------------------------
plt.subplot(4,3,7)  # , yscale="log")
u *= unit_length_in_si ** 2 / unit_time_in_si ** 2
u /= a ** (3 * (gas_gamma - 1.0))
T = (gas_gamma - 1.0) * u * mH_in_kg / k_in_J_K
print("z = {0:.2f}, T_avg = {1:.2f}".format(redshift, T.mean()))
plt.plot(x, T, **scatter_props)
plt.plot(x_s, T_s_x_upper, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(x_s, T_s_x_lower, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Temperature}}~T~{\\rm{(K)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

## Pressure profile -------------------------
plt.subplot(4,3,8)  # , yscale="log")
P /= a ** (3 * gas_gamma)
P *= unit_rho_in_si  ## not unit pressure, as the k_B/m_H changes things
plt.plot(x, P, **scatter_props)
plt.plot(x_s, P_s_x_upper, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(x_s, P_s_x_lower, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Pressure}}~P~{\\rm{(Pa)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

## Magnetic field profile -------------------------
plt.subplot(4,3,10)
Bx *= unit_B_in_si
Bx /= a**(3*gas_gamma/2)
Bx *= (a/a_i)**2   # for comoving
plt.plot(x, Bx*1e9, **scatter_props)
plt.plot(x_s, Bx_s_x_upper*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(x_s, Bx_s_x_lower*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Comoving~magnetic~field}}~B_x~{\\rm{(nT)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

plt.subplot(4,3,11)
plt.plot(y, Bx*1e9, **scatter_props)
plt.plot(y_s, Bx_s_y_upper*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(y_s, Bx_s_y_lower*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Comoving~magnetic~field}}~B_x~{\\rm{(nT)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~y~{\\rm{(Mpc)}}$", labelpad=0)

plt.subplot(4,3,12)
plt.plot(z, Bx*1e9, **scatter_props)
plt.plot(z_s, Bx_s_z_upper*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
plt.plot(z_s, Bx_s_z_lower*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
plt.ylabel("${\\rm{Comoving~magnetic~field}}~B_x~{\\rm{(nT)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~z~{\\rm{(Mpc)}}$", labelpad=0)
##

# Information -------------------------------------
plt.subplot(4,3,9, frameon=False)

text_fontsize = 5

plt.text(
    -0.45, 0.9, "Zeldovich pancake in 3D at z=%.2f " % (redshift), fontsize=text_fontsize
)
plt.text(
    -0.45,
    0.8,
    "adiabatic index $\\gamma=%.2f$, viscosity $\\alpha=%.2f$" % (gas_gamma, alpha),
    fontsize=text_fontsize,
)
plt.text(-0.45, 0.7, "magnetic diffusivity $\\eta=%.2f$" % mag_eta, fontsize=text_fontsize)
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

plt.savefig(f"{path}ZeldovichPancake_%.4d.png" % snap, dpi=200)
