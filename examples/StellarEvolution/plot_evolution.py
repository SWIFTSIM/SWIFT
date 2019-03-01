###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 ##############################################################################

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py
import numpy as np

def find_indices(a,b):
        result = np.zeros(len(b))
        for i in range(len(b)):
                result[i] = ((np.where(a == b[i]))[0])[0]

        return result


# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.1,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.1,
'figure.subplot.top'     : 0.95,
'figure.subplot.wspace'  : 0.2,
'figure.subplot.hspace'  : 0.2,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


n_snapshots = 52
n_particles_to_plot = 500

# Read the simulation data
sim = h5py.File("stellar_evolution_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

# Units
unit_length_in_cgs = sim["/Units"].attrs["Unit length in cgs (U_L)"]
unit_mass_in_cgs = sim["/Units"].attrs["Unit mass in cgs (U_M)"]
unit_time_in_cgs = sim["/Units"].attrs["Unit time in cgs (U_t)"]
unit_temp_in_cgs = sim["/Units"].attrs["Unit temperature in cgs (U_T)"]
unit_vel_in_cgs = unit_length_in_cgs / unit_time_in_cgs
unit_energy_in_cgs = unit_mass_in_cgs * unit_vel_in_cgs * unit_vel_in_cgs
unit_length_in_si = 0.01 * unit_length_in_cgs
unit_mass_in_si = 0.001 * unit_mass_in_cgs
unit_time_in_si = unit_time_in_cgs
unit_density_in_cgs = unit_mass_in_cgs*unit_length_in_cgs**-3
unit_pressure_in_cgs = unit_mass_in_cgs/unit_length_in_cgs*unit_time_in_cgs**-2
unit_int_energy_in_cgs = unit_energy_in_cgs/unit_mass_in_cgs
unit_entropy_in_cgs = unit_energy_in_cgs/unit_temp_in_cgs
Myr_in_cgs = 3.154e13
Msun_in_cgs = 1.989e33

pos = sim["/PartType0/Coordinates"][:,:]
x = pos[:,0] - boxSize / 2
y = pos[:,1] - boxSize / 2
z = pos[:,2] - boxSize / 2
vel = sim["/PartType0/Velocities"][:,:]
r = sqrt(x**2 + y**2 + z**2)
v_r = (x * vel[:,0] + y * vel[:,1] + z * vel[:,2]) / r
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]
mass = sim["/PartType0/Masses"][:]
IDs = sim["/PartType0/ParticleIDs"][:]

index = argsort(r)
part_IDs_to_plot = zeros(n_particles_to_plot)
part_IDs_to_plot = np.sort(IDs[index[0:n_particles_to_plot]])

masses_to_plot = zeros((n_particles_to_plot, n_snapshots))
v_r_to_plot = zeros((n_particles_to_plot, n_snapshots))
metallicities_to_plot = zeros((n_particles_to_plot, n_snapshots))
internal_energies_to_plot = zeros((n_particles_to_plot, n_snapshots))
t = zeros(n_snapshots)

for i in range(n_snapshots):
	print("reading snapshot "+str(i))
	# Read the simulation data
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	t[i] = sim["/Header"].attrs["Time"][0]
	
	pos = sim["/PartType0/Coordinates"][:,:]
	x = pos[:,0] - boxSize / 2
	y = pos[:,1] - boxSize / 2
	z = pos[:,2] - boxSize / 2
	vel = sim["/PartType0/Velocities"][:,:]
	r = sqrt(x**2 + y**2 + z**2)
	v_r = (x * vel[:,0] + y * vel[:,1] + z * vel[:,2]) / r
	u = sim["/PartType0/InternalEnergy"][:]
	S = sim["/PartType0/Entropy"][:]
	P = sim["/PartType0/Pressure"][:]
	rho = sim["/PartType0/Density"][:]
	mass = sim["/PartType0/Masses"][:]
	metallicity = sim["/PartType0/Metallicity"][:]
	internal_energy = sim["/PartType0/InternalEnergy"][:]
	IDs = sim["/PartType0/ParticleIDs"][:]
	indices = (find_indices(IDs,part_IDs_to_plot)).astype(int)
	masses_to_plot[:,i] = mass[indices[:]]
	v_r_to_plot[:,i] = v_r[indices[:]]
	metallicities_to_plot[:,i] = metallicity[indices[:]]
	internal_energies_to_plot[:,i] = internal_energy[indices[:]]
	

# Plot the interesting quantities
figure()

# Radial velocity --------------------------------
subplot(221)
for j in range(n_particles_to_plot):
	plot(t * unit_time_in_cgs / Myr_in_cgs, v_r_to_plot[j,:] * unit_vel_in_cgs, linewidth=0.5, color='k', ms=0.5, alpha=0.1)
xlabel("${\\rm{Time}} (Myr)$", labelpad=0)
ylabel("${\\rm{Radial~velocity}}~v_r (cm \cdot s^{-1})$", labelpad=0)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Masses --------------------------------
subplot(222)
for j in range(n_particles_to_plot):
	plot(t * unit_time_in_cgs, masses_to_plot[j,:] * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', ms=0.5, alpha=0.1)
xlabel("${\\rm{Time}} (Myr)$", labelpad=0)
ylabel("${\\rm{Mass}} (Msun)$", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ylim(0.9*np.max(masses_to_plot)*unit_mass_in_cgs/Msun_in_cgs, 1.1*np.max(masses_to_plot)*unit_mass_in_cgs/Msun_in_cgs)

# Metallicities --------------------------------
subplot(223)
for j in range(n_particles_to_plot):
	plot(t * unit_time_in_cgs, metallicities_to_plot[j,:] , linewidth=0.5, color='k', ms=0.5, alpha=0.1)
xlabel("${\\rm{Time}} (Myr)$", labelpad=0)
ylabel("${\\rm{Metallicity}} ()$", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ylim(0.9*np.max(masses_to_plot)*unit_mass_in_cgs/Msun_in_cgs, 1.1*np.max(masses_to_plot)*unit_mass_in_cgs/Msun_in_cgs)

# Internal energy --------------------------------
subplot(224)
for j in range(n_particles_to_plot):
	plot(t * unit_time_in_cgs, internal_energies_to_plot[j,:] * unit_energy_in_cgs / unit_mass_in_cgs, linewidth=0.5, color='k', ms=0.5, alpha=0.1)
xlabel("${\\rm{Time}} (Myr)$", labelpad=0)
ylabel("${\\rm{Internal energy}} (erg \cdot g^{-1})$", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ylim(0.9*np.max(masses_to_plot)*unit_mass_in_cgs/Msun_in_cgs, 1.1*np.max(masses_to_plot)*unit_mass_in_cgs/Msun_in_cgs)

savefig("time_evolution.png", dpi=200)




