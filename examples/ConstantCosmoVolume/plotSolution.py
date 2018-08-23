################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
T_i = 100.           # Initial temperature of the gas (in K)
z_c = 1.             # Redshift of caustic formation (non-linear collapse)
z_i = 100.           # Initial redshift
gas_gamma = 5./3.    # Gas adiabatic index

# Physical constants needed for internal energy to temperature conversion
kB_in_SI = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py
import os.path

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.06,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.06,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.21,
'figure.subplot.hspace'  : 0.13,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})

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
git = sim["Code"].attrs["Git Revision"]
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
m_gas = sim["/PartType0/Masses"][0]
N = sim["/Header"].attrs["NumPart_Total"][0]
sim.close()

# Expected comoving quantities
rho_0 = N * m_gas / boxSize**3
u_0 = kB_in_SI * T_i / (gas_gamma - 1.) / mH_in_kg
u_0 *= 1e-6 # conversion to internal units
u_0 *= a**(-3*(1-gas_gamma))
S_0 = (gas_gamma - 1.) * u_0 * rho_0**(-(gas_gamma - 1.))

# Mean quantities over time
z = np.zeros(119)
a = np.zeros(119)
S_mean = np.zeros(119)
S_std = np.zeros(119)
u_mean = np.zeros(119)
u_std = np.zeros(119)
P_mean = np.zeros(119)
P_std = np.zeros(119)
rho_mean = np.zeros(119)
rho_std = np.zeros(119)

vx_mean = np.zeros(119)
vy_mean = np.zeros(119)
vz_mean = np.zeros(119)
vx_std = np.zeros(119)
vy_std = np.zeros(119)
vz_std = np.zeros(119)

for i in range(119):
    sim = h5py.File("box_%04d.hdf5"%i, "r")

    z[i] = sim["/Cosmology"].attrs["Redshift"][0]
    a[i] = sim["/Cosmology"].attrs["Scale-factor"][0]
    
    S = sim["/PartType0/Entropy"][:]
    S_mean[i] = np.mean(S)
    S_std[i] = np.std(S)
    
    u = sim["/PartType0/InternalEnergy"][:]
    u_mean[i] = np.mean(u)
    u_std[i] = np.std(u)

    P = sim["/PartType0/Pressure"][:]
    P_mean[i] = np.mean(P)
    P_std[i] = np.std(P)

    rho = sim["/PartType0/Density"][:]
    rho_mean[i] = np.mean(rho)
    rho_std[i] = np.std(rho)

    v = sim["/PartType0/Velocities"][:,:]
    vx_mean[i] = np.mean(v[:,0])
    vy_mean[i] = np.mean(v[:,1])
    vz_mean[i] = np.mean(v[:,2])
    vx_std[i] = np.std(v[:,0])
    vy_std[i] = np.std(v[:,1])
    vz_std[i] = np.std(v[:,2])
    
# Move to physical quantities
rho_mean_phys = rho_mean / a**3
u_mean_phys = u_mean / a**(3*(gas_gamma - 1.))
S_mean_phys = S_mean

# Solution in physical coordinates
#T_solution = np.ones(T) / a

figure()

# Density evolution --------------------------------
subplot(231)#, yscale="log")
semilogx(a, rho_mean / rho_0, '-', color='r', lw=1)
semilogx([1e-10, 1e10], np.ones(2), '-', color='0.6', lw=1.)
semilogx([1e-10, 1e10], np.ones(2)*0.99, '--', color='0.6', lw=1.)
semilogx([1e-10, 1e10], np.ones(2)*1.01, '--', color='0.6', lw=1.)
text(1e-2, 1.0105, "+1\\%", color='0.6', fontsize=9)
text(1e-2, 0.9895, "-1\\%", color='0.6', fontsize=9, va="top")
text(1e-2, 1.015, "$\\rho_0=%.3f$"%rho_0) 
ylim(0.98, 1.02)
xlim(8e-3, 1.1)
xlabel("${\\rm Scale-factor}$", labelpad=0.)
ylabel("${\\rm Comoving~density}~\\rho / \\rho_0$", labelpad=0.)

# Thermal energy evolution --------------------------------
subplot(232)#, yscale="log")
semilogx(a, u_mean / u_0, '-', color='r', lw=1)
semilogx([1e-10, 1e10], np.ones(2), '-', color='0.6', lw=1.)
semilogx([1e-10, 1e10], np.ones(2)*0.99, '--', color='0.6', lw=1.)
semilogx([1e-10, 1e10], np.ones(2)*1.01, '--', color='0.6', lw=1.)
text(1e-2, 1.0105, "+1\\%", color='0.6', fontsize=9)
text(1e-2, 0.9895, "-1\\%", color='0.6', fontsize=9, va="top")
text(1e-2, 1.015, "$u_0=%.3e$"%(u_0)) 
ylim(0.98, 1.02)
xlim(8e-3, 1.1)
xlabel("${\\rm Scale-factor}$", labelpad=0.)
ylabel("${\\rm Comoving~internal~energy}~u / u_0$", labelpad=0.)

# Entropy evolution --------------------------------
subplot(233)#, yscale="log")
semilogx(a, S_mean / S_0, '-', color='r', lw=1)
semilogx([1e-10, 1e10], np.ones(2), '-', color='0.6', lw=1.)
semilogx([1e-10, 1e10], np.ones(2)*0.99, '--', color='0.6', lw=1.)
semilogx([1e-10, 1e10], np.ones(2)*1.01, '--', color='0.6', lw=1.)
text(1e-2, 1.0105, "+1\\%", color='0.6', fontsize=9)
text(1e-2, 0.9895, "-1\\%", color='0.6', fontsize=9, va="top")
text(1e-2, 1.015, "$A_0=%.3e$"%(S_0)) 
ylim(0.98, 1.02)
xlim(8e-3, 1.1)
xlabel("${\\rm Scale-factor}$", labelpad=0.)
ylabel("${\\rm Comoving~entropy}~A / A_0$", labelpad=0.)

# Peculiar velocity evolution ---------------------
subplot(234)
semilogx(a, vx_mean, '-', color='r', lw=1)
semilogx(a, vy_mean, '-', color='g', lw=1)
semilogx(a, vz_mean, '-', color='b', lw=1)
xlabel("${\\rm Scale-factor}$", labelpad=0.)
ylabel("${\\rm Peculiar~velocity~mean}$", labelpad=-5.)

# Peculiar velocity evolution ---------------------
subplot(235)
semilogx(a, vx_std, '--', color='r', lw=1)
semilogx(a, vy_std, '--', color='g', lw=1)
semilogx(a, vz_std, '--', color='b', lw=1)
xlabel("${\\rm Scale-factor}$", labelpad=0.)
ylabel("${\\rm Peculiar~velocity~std-dev}$", labelpad=0.)


# Information -------------------------------------
subplot(236, frameon=False)

plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

savefig("ConstantBox_comoving.png", dpi=200)



figure()

# Density evolution --------------------------------
subplot(231)#, yscale="log")
loglog(a, rho_mean_phys, '-', color='r', lw=1)
xlabel("${\\rm Scale-factor}$")
ylabel("${\\rm Physical~density}$")

# Thermal energy evolution --------------------------------
subplot(232)#, yscale="log")
loglog(a, u_mean_phys, '-', color='r', lw=1)
xlabel("${\\rm Scale-factor}$")
ylabel("${\\rm Physical~internal~energy}$")

# Entropy evolution --------------------------------
subplot(233)#, yscale="log")
semilogx(a, S_mean_phys, '-', color='r', lw=1)
xlabel("${\\rm Scale-factor}$")
ylabel("${\\rm Physical~entropy}$")

# Information -------------------------------------
subplot(236, frameon=False)

plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

savefig("ConstantBox_physical.png", dpi=200)


