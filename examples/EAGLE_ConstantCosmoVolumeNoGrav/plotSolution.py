################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2019 Stefan Arridge (stefan.arridge@durham.ac.uk)
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
T_i = 1000.           # Initial temperature of the gas (in K)
#z_c = 1.             # Redshift of caustic formation (non-linear collapse)
z_i = 100.           # Initial redshift
gas_gamma = 5./3.    # Gas adiabatic index

# Physical constants needed for internal energy to temperature conversion
kB_in_SI = 1.38064852e-23
mH_in_kg = 1.6737236e-27

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
sim = h5py.File("eagle_cooling_box_0000.hdf5", "r")
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

# Mean quantities over time
z = np.zeros(119)
a = np.zeros(119)
T_mean = np.zeros(119)
T_std = np.zeros(119)

for i in range(119):
    sim = h5py.File("eagle_cooling_box_%04d.hdf5"%i, "r")

    z[i] = sim["/Cosmology"].attrs["Redshift"][0]
    a[i] = sim["/Cosmology"].attrs["Scale-factor"][0]
    
    T = sim["/PartType0/Temperature"][0]/1.0e4
    T_mean[i] = np.mean(T)
    T_std[i] = np.std(T)

# Read in observational data

data_schaye = np.genfromtxt("./datasets/schaye_et_al_2000_thermal_history.dat",skip_header = 2)
data_walther = np.genfromtxt("./datasets/walther_et_al_2019_thermal_history.dat",skip_header = 2)

data_schaye = data_schaye.T
data_walther = data_walther.T
print(np.shape(data_schaye))
print(data_schaye[0])

schaye_z_lower_error = data_schaye[0] - data_schaye[1]
schaye_z_upper_error = data_schaye[2] - data_schaye[0]
schaye_T_lower_error = data_schaye[3] - data_schaye[4]
schaye_T_upper_error = data_schaye[5] - data_schaye[3]
walther_T_lower_error = data_walther[1] - data_walther[2]
walther_T_upper_error = data_walther[3] - data_walther[1]

# Make plot of temperature evolution  --------------------------------
fig, ax = plt.subplots()
ax.fill_between(z,T_mean - T_std,T_mean + T_std,alpha = 0.5)
ax.plot(z,T_mean,label = "Simulation")
ax.errorbar(data_schaye[0],data_schaye[3], xerr = [schaye_z_lower_error,schaye_z_upper_error],yerr = [schaye_T_lower_error,schaye_T_upper_error], fmt = 'ko', label = "Schaye+ 2000",zorder = 20,capsize = 4.0,capthick = 1.0,alpha = 0.9)
ax.errorbar(data_walther[0],data_walther[1],yerr = [walther_T_lower_error,walther_T_upper_error], fmt = 'rd', label = "Walther+ 2019",zorder = 30,capsize = 4.0,capthick = 1.0,alpha = 0.7)
ax.set_xlim(0.0,10.0)
ax.set_ylim(0.0,3.0)
ax.set_xlabel("z")
ax.set_ylabel(r"$T_0\,/\,10^4\,\,K$")
ax.legend(loc = 0)
fig.tight_layout()
fig.savefig("thermal_history.pdf",format = "pdf")


