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
 'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.145,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.11,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


import numpy as np
import h5py as h5
import sys

# File containing the total energy
stats_filename = "./energy.txt"

# First snapshot
snap_filename = "coolingBox_0000.hdf5"

# Some constants in cgs units
k_b = 1.38E-16 #boltzmann
m_p = 1.67e-24 #proton mass

# Initial conditions set in makeIC.py
T_init = 1.0e5

# Read the initial state of the gas
f = h5.File(snap_filename,'r')
rho = np.mean(f["/PartType0/Density"])
pressure = np.mean(f["/PartType0/Pressure"])

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]

# Read the properties of the cooling function
parameters = f["Parameters"]
cooling_lambda = float(parameters.attrs["LambdaCooling:lambda_cgs"])
min_T = float(parameters.attrs["LambdaCooling:minimum_temperature"])
mu = float(parameters.attrs["LambdaCooling:mean_molecular_weight"])
X_H = float(parameters.attrs["LambdaCooling:hydrogen_mass_abundance"])

# Read the adiabatic index
gamma = float(f["HydroScheme"].attrs["Adiabatic index"])

print "Initial density :", rho
print "Initial pressure:", pressure
print "Adiabatic index :", gamma

# Read energy and time arrays
array = np.genfromtxt(stats_filename,skip_header = 1)
time = array[:,0]
total_mass = array[:,1]
total_energy = array[:,2]
kinetic_energy = array[:,3]
internal_energy = array[:,4]
radiated_energy = array[:,8]
initial_energy = total_energy[0]

# Conversions to cgs
rho_cgs = rho * unit_mass / (unit_length)**3
time_cgs = time * unit_time
total_energy_cgs = total_energy / total_mass[0] * unit_length**2 / (unit_time)**2
kinetic_energy_cgs = kinetic_energy / total_mass[0] * unit_length**2 / (unit_time)**2
internal_energy_cgs = internal_energy / total_mass[0] * unit_length**2 / (unit_time)**2
radiated_energy_cgs = radiated_energy / total_mass[0] * unit_length**2 / (unit_time)**2  

# Find the energy floor
u_floor_cgs = k_b * min_T / (mu * m_p * (gamma - 1.))

# Find analytic solution
initial_energy_cgs = initial_energy/total_mass[0] * unit_length**2 / (unit_time)**2 
n_H_cgs = X_H * rho_cgs / m_p
du_dt_cgs = -cooling_lambda * n_H_cgs**2 / rho_cgs
cooling_time_cgs = (initial_energy_cgs/(-du_dt_cgs))[0]
analytic_time_cgs = np.linspace(0, cooling_time_cgs * 1.8, 1000)
u_analytic_cgs = du_dt_cgs*analytic_time_cgs + initial_energy_cgs
u_analytic_cgs[u_analytic_cgs < u_floor_cgs] = u_floor_cgs

print "Cooling time:", cooling_time_cgs, "[s]"

# Read snapshots
u_snapshots_cgs = zeros(25)
t_snapshots_cgs = zeros(25)
for i in range(25):
    snap = h5.File("coolingBox_%0.4d.hdf5"%i,'r')
    u_snapshots_cgs[i] = sum(snap["/PartType0/InternalEnergy"][:] * snap["/PartType0/Masses"][:])  / total_mass[0] * unit_length**2 / (unit_time)**2
    t_snapshots_cgs[i] = snap["/Header"].attrs["Time"] * unit_time


figure()
plot(time_cgs, total_energy_cgs, 'r-', lw=1.6, label="Gas total energy")
plot(t_snapshots_cgs, u_snapshots_cgs, 'rD', ms=3)
plot(time_cgs, radiated_energy_cgs, 'g-', lw=1.6, label="Radiated energy")
plot(time_cgs, total_energy_cgs + radiated_energy_cgs, 'b-', lw=0.6, label="Gas total + radiated")

plot(analytic_time_cgs, u_analytic_cgs, '--', color='k', alpha=0.8, lw=1.0, label="Analytic solution")

legend(loc="upper right", fontsize=8, frameon=False, handlelength=3, ncol=1)
xlabel("${\\rm{Time~[s]}}$", labelpad=0)
ylabel("${\\rm{Energy~[erg]}}$")
xlim(0, 1.5*cooling_time_cgs)
ylim(0, 1.5*u_analytic_cgs[0])

savefig("energy.png", dpi=200)


