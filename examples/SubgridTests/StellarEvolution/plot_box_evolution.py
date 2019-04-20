###############################################################################
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
 ##############################################################################

# Script used to plot time evolution of gas particle properties. Intended to 
# compare result of feedback due to one star placed in centre of uniform box
# of gas with output from EAGLE feedback test. Can also use as input output
# from SWIFT feedback test (tests/testFeedback) with the appropriate change
# to filepath.

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py
import numpy as np
import glob
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
'figure.subplot.left'    : 0.05,
'figure.subplot.right'   : 0.995,
'figure.subplot.bottom'  : 0.06,
'figure.subplot.top'     : 0.92,
'figure.subplot.wspace'  : 0.25,
'figure.subplot.hspace'  : 0.2,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


# Number of snapshots and elements
newest_snap_name = max(glob.glob('stellar_evolution_*.hdf5'))#, key=os.path.getctime)
n_snapshots = int(newest_snap_name.replace('stellar_evolution_','').replace('.hdf5','')) + 1
n_elements = 9

# Read the simulation data
sim = h5py.File("stellar_evolution_0000.hdf5", "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]
stellar_mass = sim["/PartType4/Masses"][0]
E_SNIa_cgs = double(sim["/Parameters"].attrs["EAGLEFeedback:SNIa_energy_erg"])

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
Gyr_in_cgs = 1e9 * 365. * 24 * 3600.
Msun_in_cgs = 1.98848e33

# Declare arrays to store SWIFT data
swift_box_gas_mass = zeros(n_snapshots)
swift_box_star_mass = zeros(n_snapshots)
swift_box_gas_metal_mass = zeros(n_snapshots)
swift_element_mass = zeros((n_snapshots,n_elements))
t = zeros(n_snapshots)

# Read data from snapshots
for i in range(n_snapshots):
	#print("reading snapshot "+str(i))
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	t[i] = sim["/Header"].attrs["Time"][0]
	
	masses = sim["/PartType0/Masses"][:]
	swift_box_gas_mass[i] = np.sum(masses)

        Z_star = sim["/PartType4/Metallicity"][0]
	star_masses = sim["/PartType4/Masses"][:]
	swift_box_star_mass[i] = np.sum(star_masses)

	metallicities = sim["/PartType0/Metallicity"][:]
	swift_box_gas_metal_mass[i] = np.sum(metallicities * masses)

	element_abundances = sim["/PartType0/ElementAbundance"][:][:]
	for j in range(n_elements):
		swift_element_mass[i,j] = np.sum(element_abundances[:,j] * masses)

        sim.close()

# Read expected yields from EAGLE. Choose which file to use based on metallicity used when
# running SWIFT (can be specified in yml file)
filename = "./StellarEvolutionSolution/Z_%.4f/StellarEvolutionTotal.txt"%Z_star

# Read EAGLE test output
with open(filename) as f:
	eagle_categories = f.readline()
	eagle_data = f.readlines()

eagle_data = [x.strip() for x in eagle_data]

# Declare arrays to store EAGLE test output
eagle_time_Gyr = zeros(len(eagle_data))
eagle_total_mass = zeros(len(eagle_data))
eagle_total_metal_mass = zeros(len(eagle_data))
eagle_total_element_mass = zeros((len(eagle_data),n_elements))

# Populate arrays with data from EAGLE test output
i = 0
for line in eagle_data:
	enrich_to_date = line.split(' ')
	eagle_time_Gyr[i] = float(enrich_to_date[0])
	eagle_total_mass[i] = float(enrich_to_date[1]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
	eagle_total_metal_mass[i] = float(enrich_to_date[2]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs 
	for j in range(n_elements):
		eagle_total_element_mass[i,j] = float(enrich_to_date[3+j]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
	i += 1

# Plot the interesting quantities
figure()

suptitle("Star metallicity Z = %.4f"%Z_star)

# Box gas mass --------------------------------
subplot(221)
plot(t[1:] * unit_time_in_cgs / Gyr_in_cgs, (swift_box_gas_mass[1:] - swift_box_gas_mass[0])* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', marker = "*", ms=0.5, label='swift')
plot(eagle_time_Gyr[1:],eagle_total_mass[:-1],linewidth=0.5,color='r',label='eagle test total', ls='--')
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total gas particle mass (Msun)", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
legend()

# Box star mass --------------------------------
subplot(222)
plot(t * unit_time_in_cgs / Gyr_in_cgs, (swift_box_star_mass)* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', marker = "*", ms=0.5, label='swift')
plot(eagle_time_Gyr[1:], swift_box_star_mass[0] * unit_mass_in_cgs / Msun_in_cgs - eagle_total_mass[:-1],linewidth=0.5,color='r',label='eagle test total')
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total star particle mass (Msun)", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
legend()

# Box gas element  mass --------------------------------
colours = ['k','r','g','b','c','y','m','skyblue','plum']
element_names = ['H','He','C','N','O','Ne','Mg','Si','Fe']
subplot(223)
for j in range(n_elements):
	plot(t[1:] * unit_time_in_cgs / Gyr_in_cgs, (swift_element_mass[1:,j] - swift_element_mass[0,j]) * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color=colours[j], ms=0.5, label=element_names[j])
	plot(eagle_time_Gyr[1:],eagle_total_element_mass[:-1,j],linewidth=1,color=colours[j],linestyle='--')
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in element mass of gas particles (Msun)", labelpad=2)
xscale("log")
yscale("log")
legend(bbox_to_anchor=(1.005, 1.), ncol=1, fontsize=8, handlelength=1)

# Box gas metal mass --------------------------------
subplot(224)
plot(t[1:] * unit_time_in_cgs / Gyr_in_cgs, (swift_box_gas_metal_mass[1:] - swift_box_gas_metal_mass[0])* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', marker = "*", ms=0.5, label='swift')
plot(eagle_time_Gyr[1:],eagle_total_metal_mass[:-1],linewidth=0.5,color='r',label='eagle test')
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total metal mass of gas particles (Msun)", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))

savefig("box_evolution_Z_%.4f.png"%(Z_star), dpi=200)




