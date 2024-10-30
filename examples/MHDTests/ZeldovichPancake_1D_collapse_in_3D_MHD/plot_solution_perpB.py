### This file is part of SWIFT.

# Computes the analytical solution of the Zeldovich pancake and compares with
# the simulation result, including a perpendicular magnetic field

# Parameters
B_0 = 1.0e-9 # Initial magnetic field strength (in T)
T_i = 100.0  # Initial temperature of the gas (in K)
z_c = 1.0  # Redshift of caustic formation (non-linear collapse)
z_i = 100.0  # Initial redshift
a_c = 1/(1+z_c)
a_i = 1/(1+z_i)
numPart_1D = 32  # Number of particles along each dimension

# Physical constants needed for internal energy to temperature conversion
k_in_J_K = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys
import os.path

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])
# You can input a different directory where a snapshot is saved; if you don't, use the current one
try:
    path = str(sys.argv[2])
except:
    path = ''
filepath = f"{path}zeldovichPancake_%04d.hdf5" % snap

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
v = sim["/PartType0/Velocities"][:, 0]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]
m = sim["/PartType0/Masses"][:]
x -= 0.5 * boxSize

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

# Solution taken from Springel 2010. Eqs. 127 - 130
q = np.linspace(-0.5 * lambda_i, 0.5 * lambda_i, 256)
k_i = 2.0 * np.pi / lambda_i
zfac = (1.0 + z_c) / (1.0 + redshift)
x_s = q - zfac * np.sin(k_i * q) / k_i
rho_s = rho_0 / (1.0 - zfac * np.cos(k_i * q))
v_s = -H_0 * (1.0 + z_c) / np.sqrt(1.0 + redshift) * np.sin(k_i * q) / k_i
T_s = T_i * (((1.0 + redshift) / (1.0 + z_i)) ** 3.0 * rho_s / rho_0) ** (2.0 / 3.0)
P_s = T_s * rho_s*unit_rho_in_si * k_in_J_K / mH_in_kg / (a/a_i)**3  ## rho_s is not SI units! So multiply unit

# Magnetic field solution
rho_initial = rho_0 / (1.0 - (1.0 + z_c) / (1.0 + z_i) * np.cos(k_i * q))
By_s = B_0 * rho_s / rho_initial
Bx_s = np.zeros_like(rho_s)
Bz_s = np.zeros_like(rho_s)

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

# Velocity profile --------------------------------
plt.subplot(331)
plt.plot(x, v, **scatter_props)
if redshift > z_c:
    plt.plot(x_s, v_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)
plt.ylabel("${\\rm{Peculiar~velocity}}~v_x~{\\rm{(km/s)}}$", labelpad=0)

# Density profile --------------------------------
plt.subplot(332)  # , yscale="log")
if redshift > z_c:
    plt.plot(x, rho / rho_0, **scatter_props)
    plt.plot(x_s, rho_s / rho_0, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.ylabel("${\\rm{Comoving~density}}~\\rho / \\rho_0$", labelpad=0)
else:
    plt.plot(x, np.log10(rho / rho_0), **scatter_props)
    plt.ylabel("${\\rm{Comoving~density}}~\\log_{10}(\\rho / \\rho_0)$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

# Temperature profile -------------------------
plt.subplot(334)  # , yscale="log")
u *= unit_length_in_si ** 2 / unit_time_in_si ** 2
u /= a ** (3 * (gas_gamma - 1.0))
T = (gas_gamma - 1.0) * u * mH_in_kg / k_in_J_K
print("z = {0:.2f}, T_avg = {1:.2f}".format(redshift, T.mean()))
if redshift > z_c:
    plt.plot(x, T, **scatter_props)
    plt.plot(x_s, T_s, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.ylabel("${\\rm{Temperature}}~T~{\\rm{(K)}}$", labelpad=0)
else:
    plt.plot(x, np.log10(T), **scatter_props)
    plt.ylabel("${\\rm{Temperature}}~\\log_{10}(T)~{\\rm{(K)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

# Pressure profile -------------------------
plt.subplot(335)  # , yscale="log")
P /= a ** (3 * gas_gamma)
P *= unit_rho_in_si  ## not unit pressure, as the k_B/m_H changes things
if redshift > z_c:
    plt.plot(x, P, **scatter_props)
    plt.plot(x_s, P_s, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.ylabel("${\\rm{Pressure}}~P~{\\rm{(Pa)}}$", labelpad=0)
else:
    plt.plot(x, np.log10(P), **scatter_props)
    plt.ylabel("${\\rm{Pressure}}~\\log_{10}(P)~{\\rm{(Pa)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

# Magnetic field profile -------------------------
plt.subplot(337)
Bx *= unit_B_in_si
Bx /= a**(3*gas_gamma/2)   # In ODI
## In FDI, need to add: Bx *= a**0.5
## In VP,  need to add: Bx *= a**0.5  AND  Bx /= a_i
Bx *= (a/a_i)**2   # for comoving

if redshift > z_c:
    plt.plot(x, Bx*1e9, **scatter_props)
    plt.plot(x_s, Bx_s*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.ylabel("${\\rm{Comoving~magnetic~field}}~B_x~{\\rm{(nT)}}$", labelpad=0)
else:
    plt.plot(x, np.log10(Bx*1e9), **scatter_props)
    plt.ylabel("${\\rm{Comoving~magnetic~field}}$ \n $\\log_{10}(B_x)~{\\rm{(nT)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

plt.subplot(338)
By *= unit_B_in_si
By /= a**(3*gas_gamma/2)   # In ODI
## In FDI, need to add: By *= a**0.5
## In VP,  need to add: By *= a**0.5  AND  By /= a_i
By *= (a/a_i)**2   # for comoving

if redshift > z_c:
    plt.plot(x, By*1e9, **scatter_props)
    plt.plot(x_s, By_s*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.ylabel("${\\rm{Comoving~magnetic~field}}~B_y~{\\rm{(nT)}}$", labelpad=0)
else:
    plt.plot(x, np.log10(By*1e9), **scatter_props)
    plt.ylabel("${\\rm{Comoving~magnetic~field}}$ \n $\\log_{10}(B_y)~{\\rm{(nT)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)

plt.subplot(339)
Bz *= unit_B_in_si
Bz /= a**(3*gas_gamma/2)   # In ODI
## In FDI, need to add: Bz *= a**0.5
## In VP,  need to add: Bz *= a**0.5  AND  Bz /= a_i
Bz *= (a/a_i)**2   # for comoving

if redshift > z_c:
    plt.plot(x, Bz*1e9, **scatter_props)
    plt.plot(x_s, Bz_s*1e9, "--", color=line_color, alpha=0.8, lw=1.2)
    plt.ylabel("${\\rm{Comoving~magnetic~field}}~B_z~{\\rm{(nT)}}$", labelpad=0)
else:
    plt.plot(x, np.log10(Bz*1e9), **scatter_props)
    plt.ylabel("${\\rm{Comoving~magnetic~field}}$ \n $\\log_{10}(B_z)~{\\rm{(nT)}}$", labelpad=0)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)
##

## Divergence error -------------------------
plt.subplot(336)

# Read in divB, then also apply the same scalings as for B
divB = sim["/PartType0/MagneticDivergences"]
divB *= unit_B_in_si
divB /= a**(3*gas_gamma/2)   # In ODI
divB *= (a/a_i)**2 # for comoving
## In FDI, need to add: divB *= a**0.5
## In VP,  need to add: divB *= a**0.5  AND  divB /= a_i
B = np.sqrt(Bx**2 + By**2 + Bz**2)
# Smoothing length h
h = sim["/PartType0/SmoothingLengths"]

divB_measure = np.zeros_like(divB)
if By_s[0] != 0:   # if run without field keep divB zero
    divB_measure = h * np.abs(divB) / B

plt.plot(x, divB_measure, **scatter_props)
plt.xlabel("${\\rm{Comoving~position}}~x~{\\rm{(Mpc)}}$", labelpad=0)
plt.ylabel("${h~|\\nabla \\cdot \\mathbf{B}|~/~|\\mathbf{B}|}$", labelpad=0)

# Information -------------------------------------
plt.subplot(333, frameon=False)

text_fontsize = 5

plt.text(
    -0.45, 0.9, "Zeldovich pancake at z=%.2f " % (redshift), fontsize=text_fontsize
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
