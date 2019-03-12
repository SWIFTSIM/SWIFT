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

# Computes the temperature evolution of the gas in a cosmological box

# Physical constants needed for internal energy to temperature conversion
k_in_J_K = 1.38064852e-23
mH_in_kg = 1.6737236e-27

# Number of snapshots generated
n_snapshots = 200

import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py
import os.path

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 9,
'legend.fontsize': 9,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.14,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.12,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 2.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})

# Read the simulation data
sim = h5py.File("snap_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"][0]
kernel = sim["/HydroScheme"].attrs["Kernel function"][0]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = sim["/HydroScheme"].attrs["Kernel eta"][0]
alpha = sim["/HydroScheme"].attrs["Alpha viscosity"][0]
H_mass_fraction = sim["/HydroScheme"].attrs["Hydrogen mass fraction"][0]
H_transition_temp = sim["/HydroScheme"].attrs["Hydrogen ionization transition temperature"][0]
T_initial = sim["/HydroScheme"].attrs["Initial temperature"][0]
T_minimal = sim["/HydroScheme"].attrs["Minimal temperature"][0]
git = sim["Code"].attrs["Git Revision"]
cooling_model = sim["/SubgridScheme"].attrs["Cooling Model"]

if cooling_model == "Constant Lambda":
    Lambda = sim["/SubgridScheme"].attrs["Lambda/n_H^2 [cgs]"][0]   
    
# Cosmological parameters
H_0 = sim["/Cosmology"].attrs["H0 [internal units]"][0]
gas_gamma = sim["/HydroScheme"].attrs["Adiabatic index"][0]

unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]

unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs

# Primoridal mean molecular weight as a function of temperature
def mu(T, H_frac=H_mass_fraction, T_trans=H_transition_temp):
    if T > T_trans:
        return 4. / (8. - 5. * (1. - H_frac))
    else:
        return 4. / (1. + 3. * H_frac)
    
# Temperature of some primoridal gas with a given internal energy
def T(u, H_frac=H_mass_fraction, T_trans=H_transition_temp):
    T_over_mu = (gas_gamma - 1.) * u * mH_in_kg / k_in_J_K
    ret = np.ones(np.size(u)) * T_trans

    # Enough energy to be ionized?
    mask_ionized = (T_over_mu > (T_trans+1) / mu(T_trans+1, H_frac, T_trans))
    if np.sum(mask_ionized)  > 0:
        ret[mask_ionized] = T_over_mu[mask_ionized] * mu(T_trans*10, H_frac, T_trans)

    # Neutral gas?
    mask_neutral = (T_over_mu < (T_trans-1) / mu((T_trans-1), H_frac, T_trans))
    if np.sum(mask_neutral)  > 0:
        ret[mask_neutral] = T_over_mu[mask_neutral] * mu(0, H_frac, T_trans)
        
    return ret

z = np.zeros(n_snapshots)
a = np.zeros(n_snapshots)
T_mean = np.zeros(n_snapshots)
T_std = np.zeros(n_snapshots)
T_log_mean = np.zeros(n_snapshots)
T_log_std = np.zeros(n_snapshots)
T_median = np.zeros(n_snapshots)
T_min = np.zeros(n_snapshots)
T_max = np.zeros(n_snapshots)

# Loop over all the snapshots
for i in range(n_snapshots):
    sim = h5py.File("snap_%04d.hdf5"%i, "r")

    z[i] = sim["/Cosmology"].attrs["Redshift"][0]
    a[i] = sim["/Cosmology"].attrs["Scale-factor"][0]

    u = sim["/PartType0/InternalEnergy"][:]

    # Compute the temperature
    u *= (unit_length_in_si**2 / unit_time_in_si**2)
    u /= a[i]**(3 * (gas_gamma - 1.))
    Temp = T(u)

    # Gather statistics
    T_median[i] = np.median(Temp)
    T_mean[i] = Temp.mean()
    T_std[i] = Temp.std()
    T_log_mean[i] = np.log10(Temp).mean()
    T_log_std[i] = np.log10(Temp).std()
    T_min[i] = Temp.min()
    T_max[i] = Temp.max()

# CMB evolution
a_evol = np.logspace(-3, 0, 60)
T_cmb = (1. / a_evol)**2 * 2.72

# Plot the interesting quantities
figure()
subplot(111, xscale="log", yscale="log")

fill_between(a, T_mean-T_std, T_mean+T_std, color='C0', alpha=0.1)
plot(a, T_max, ls='-.', color='C0', lw=1., label="${\\rm max}~T$")
plot(a, T_min, ls=':', color='C0', lw=1., label="${\\rm min}~T$")
plot(a, T_mean, color='C0', label="${\\rm mean}~T$", lw=1.5)
fill_between(a, 10**(T_log_mean-T_log_std), 10**(T_log_mean+T_log_std), color='C1', alpha=0.1)
plot(a, 10**T_log_mean, color='C1', label="${\\rm mean}~{\\rm log} T$", lw=1.5)
plot(a, T_median, color='C2', label="${\\rm median}~T$", lw=1.5)

legend(loc="upper left", frameon=False, handlelength=1.5)

# Cooling model
if cooling_model == "Constant Lambda":
    text(1e-2, 6e4, "$\Lambda_{\\rm const}/n_{\\rm H}^2 = %.1f\\times10^{%d}~[\\rm{cgs}]$"%(Lambda/10.**(int(log10(Lambda))), log10(Lambda)), fontsize=7)
elif cooling_model == "EAGLE":
    text(1e-2, 6e4, "EAGLE (Wiersma et al. 2009)")
elif cooling_model == b"Grackle":
    text(1e-2, 6e4, "Grackle (Smith et al. 2016)")
else:
    text(1e-2, 6e4, "No cooling")
    
# Expected lines
plot([1e-10, 1e10], [H_transition_temp, H_transition_temp], 'k--', lw=0.5, alpha=0.7)
text(2.5e-2, H_transition_temp*1.07, "$T_{\\rm HII\\rightarrow HI}$", va="bottom", alpha=0.7, fontsize=8)
plot([1e-10, 1e10], [T_minimal, T_minimal], 'k--', lw=0.5, alpha=0.7)
text(1e-2, T_minimal*0.8, "$T_{\\rm min}$", va="top", alpha=0.7, fontsize=8)
plot(a_evol, T_cmb, 'k--', lw=0.5, alpha=0.7)
text(a_evol[20], T_cmb[20]*0.55, "$(1+z)^2\\times T_{\\rm CMB,0}$", rotation=-34, alpha=0.7, fontsize=8, va="top", bbox=dict(facecolor='w', edgecolor='none', pad=1.0, alpha=0.9))


redshift_ticks = np.array([0., 1., 2., 5., 10., 20., 50., 100.])
redshift_labels = ["$0$", "$1$", "$2$", "$5$", "$10$", "$20$", "$50$", "$100$"]
a_ticks = 1. / (redshift_ticks + 1.)

xticks(a_ticks, redshift_labels)
minorticks_off()

xlabel("${\\rm Redshift}~z$", labelpad=0)
ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=0)
xlim(9e-3, 1.1)
ylim(5, 2.5e7)

savefig("Temperature_evolution.png", dpi=200)

