################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
gas_gamma = 5./3.    # Polytropic index
lambda_i = 1.975e24  # Wavelength of the density perturbation (in h^-1 m)
rho_0 = 1.8788e-26   # Critical density of the Universe (in h^2 kg m^-3)
T_i = 100.           # Initial temperature of the gas (in K)
H_0 = 3.24077929e-18 # Assumed value of the Hubble constant (in s^-1)
z_c = 1.             # Redshift of caustic formation (non-linear collapse)
z_i = 100.           # Initial redshift

k_in_J_K = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.045,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.05,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})

snap = int(sys.argv[1])

# Read the simulation data
sim = h5py.File("zeldovichPancake_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
redshift = sim["/Header"].attrs["Redshift"][0]
a = sim["/Header"].attrs["Scale-factor"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:,0]
v = sim["/PartType0/Velocities"][:,0]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

x_s = linspace(0., lambda_i, 1000)
k_i = 2. * pi / lambda_i
zfac = (1. + z_c) / (1. + redshift)
rho_s = rho_0 / (1 - zfac * cos(k_i * x_s))
v_s = -H_0 * (1. + z_c) / sqrt(1. + redshift) * sin(k_i * x_s) / k_i
T_s = T_i * \
  (((1. + redshift) / (1. + z_i))**3 / (1. - zfac * cos(k_i * x_s)))**(2. / 3.)
P_s = zeros(1000)
u_s = zeros(1000)
s_s = zeros(1000)

unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]

unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs

x_s /= unit_length_in_si
rho_s *= (unit_length_in_si**3 / unit_mass_in_si)
v_s *= (unit_time_in_si / unit_length_in_si)

# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
plot(x, v, '.', color='r', ms=4.0)
plot(x_s, v_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)

# Density profile --------------------------------
subplot(232)
plot(x, rho, '.', color='r', ms=4.0)
plot(x_s, rho_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)

# Pressure profile --------------------------------
subplot(233)
plot(x, P, '.', color='r', ms=4.0)
plot(x_s, P_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)

# Internal energy profile -------------------------
subplot(234)
#plot(x, u, '.', color='r', ms=4.0)
#plot(x_s, u_s, '--', color='k', alpha=0.8, lw=1.2)
u *= (unit_length_in_si**2 / unit_time_in_si**2)
u /= a**(3 * (gas_gamma - 1.))
T = (gas_gamma - 1.) * u * mH_in_kg / k_in_J_K
print T.mean()
plot(x, T, '.', color='r', ms=4.0)
plot(x_s, T_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)

# Entropy profile ---------------------------------
subplot(235)
plot(x, S, '.', color='r', ms=4.0)
plot(x_s, s_s, '--', color='k', alpha=0.8, lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Zeldovich pancake with  $\\gamma=%.3f$ in 1D at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "$z={0:.2f}$".format(redshift))
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

savefig("ZeldovichPancake.png", dpi=200)
