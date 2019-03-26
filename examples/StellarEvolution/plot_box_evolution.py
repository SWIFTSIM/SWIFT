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
import glob
import os.path

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


# Number of snapshots and elements
newest_snap_name = max(glob.glob('stellar_evolution_*.hdf5'), key=os.path.getctime)
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
Gyr_in_cgs = 3.155e16
Msun_in_cgs = 1.989e33

box_mass = zeros(n_snapshots)
box_metal_mass = zeros(n_snapshots)
element_mass = zeros((n_snapshots,n_elements))
#abundance = zeros((n_parts,n_elements,n_snapshots))
t = zeros(n_snapshots)

# Read expected yields from EAGLE
#filename = "/cosma/home/dp004/dc-bori1/Eagle/data1/StellarEvolutionTotal.txt"
#filename = "/cosma/home/dp004/dc-bori1/Eagle/data1/StellarEvolutionIa.txt"
filename = "/cosma/home/dp004/dc-bori1/Eagle/data1/StellarEvolutionII.txt"
#filename = "/cosma/home/dp004/dc-bori1/Eagle/data1/StellarEvolutionAGB.txt"

with open(filename) as f:
	eagle_categories = f.readline()
	eagle_data = f.readlines()

eagle_data = [x.strip() for x in eagle_data]

i = 0
time_Gyr = zeros(len(eagle_data))
total_mass = zeros(len(eagle_data))
total_metal_mass = zeros(len(eagle_data))
total_element_mass = zeros((len(eagle_data),n_elements))

for line in eagle_data:
	enrich_to_date = line.split(' ')
	time_Gyr[i] = float(enrich_to_date[0])
	total_mass[i] = float(enrich_to_date[1]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
	total_metal_mass[i] = float(enrich_to_date[2]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs 
	total_element_mass[i,0] = float(enrich_to_date[3]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs 
	i += 1

# Read expected yields from SWIFT test
#swift_test_filename = "../../tests/test_feedback_total.txt"
##swift_test_filename = "../../tests/test_feedback_SNIa.txt"
##swift_test_filename = "../../tests/test_feedback_SNII.txt"
##swift_test_filename = "../../tests/test_feedback_AGB.txt"
#
#with open(swift_test_filename) as f:
#	swift_categories = f.readline()
#	swift_data = f.readlines()
#
#swift_data = [x.strip() for x in swift_data]
#
#i = 0
#swift_test_time_Gyr = zeros(len(swift_data))
#swift_test_total_mass = zeros(len(swift_data))
#swift_test_total_metal_mass = zeros(len(swift_data))
#swift_test_total_element_mass = zeros((len(swift_data),n_elements))
#
#for line in swift_data:
#	enrich_to_date = line.split(' ')
#	swift_test_time_Gyr[i] = float(enrich_to_date[0])
#	swift_test_total_mass[i] = float(enrich_to_date[1]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs
#	swift_test_total_metal_mass[i] = float(enrich_to_date[2]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs 
#	swift_test_total_element_mass[i,0] = float(enrich_to_date[3]) * stellar_mass / Msun_in_cgs * unit_mass_in_cgs 
#	i += 1

# Read data from snapshots
for i in range(n_snapshots):
	print("reading snapshot "+str(i))
	# Read the simulation data
	sim = h5py.File("stellar_evolution_%04d.hdf5"%i, "r")
	t[i] = sim["/Header"].attrs["Time"][0]
	#ids = sim["/PartType0/ParticleIDs"][:]
	
	masses = sim["/PartType0/Masses"][:]
	box_mass[i] = np.sum(masses)
	#masses_sort = np.concatenate(ids,masses,axis = 0)
	#masses_sort = np.sort(masses_sort,axis = 0)
	#
	#metallicities = sim["/PartType0/Metallicity"][:]
	#box_metal_mass[i] = np.sum(metallicities * masses)
	#
	#abundance[:,:,i] = np.concatetenate((ids,sim["/PartType0/ElementAbundance"][:,:]),axis = 0)
	#abundance[:,:,i] = np.sort(abundance[:,:,i],axis = 0)
	#for j in range(n_parts):
	#	if (i > 0 && (abundance[j,0,i] - abundance[j,0,i-1]) != 0):
	#		element_mass[i,0] += (abundance[j,0,i] - abundance[j,0,i-1])
	#print(np.count_nonzero(abundance))
	#print(np.count_nonzero(element_mass[i,0]))

# Read data from output
#swift_output_file = "swift-feedback.out"
#swift_time = [0]
#swift_enrich_mass = [0]
#with open(swift_output_file) as f:
#	for line in f:
#		str_pos = line.find('stars_evolve_spart')
#		if str_pos != -1:
#			line = line[(str_pos+19):]
#			line = line.strip()
#			swift_output_data = line.split(' ')
#			swift_time.append(float(swift_output_data[0]))
#			swift_enrich_mass.append(swift_enrich_mass[-1] + float(swift_output_data[1]))

# Plot the interesting quantities
figure()

# Box mass --------------------------------
#subplot(221)
subplot(111)
plot(t[1:] * unit_time_in_cgs / Gyr_in_cgs, (box_mass[1:] - box_mass[0])* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', marker = "*", ms=0.5, label='swift')
plot(time_Gyr[1:],total_mass[:-1],linewidth=0.5,color='r',label='eagle test total')
#plot(swift_test_time_Gyr[1:],swift_test_total_mass[:-1],linewidth=0.5,color='c',label='swift test total')
plot(time_Gyr[1:],total_mass[:-1] - total_metal_mass[:-1],linewidth=0.5,color='g',label='eagle test total - metals')
#plot(swift_time,swift_enrich_mass * unit_mass_in_cgs / Msun_in_cgs,linewidth=0.5,color='b',label='swift output')
xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
ylabel("Change in total gas particle mass (Msun)", labelpad=2)
#ylabel("Ratio swift/eagle change in mass", labelpad=2)
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
legend()

# Box mass ratio --------------------------------
#subplot(122)
#plot(t[1:] * unit_time_in_cgs / Gyr_in_cgs, ((box_mass[1:] - box_mass[0])* unit_mass_in_cgs / Msun_in_cgs)/total_mass[0:64], linewidth=0.5, color='k', marker = "*", ms=0.5, label='swift')
#xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
#ylabel("Ratio swift/eagle change in mass", labelpad=2)
#ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Box metal mass --------------------------------
#subplot(222)
#plot(t * unit_time_in_cgs / Gyr_in_cgs, (box_metal_mass - box_metal_mass[0])* unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', ms=0.5, label='swift')
#plot(time_Gyr,total_metal_mass,linewidth=0.5,color='r',label='eagle test')
#xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
#ylabel("Change in total metal mass of gas particles (Msun)", labelpad=2)
#ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#legend()
#
## Box element mass --------------------------------
#subplot(223)
#plot(t * unit_time_in_cgs / Gyr_in_cgs, (element_mass[:,0] - element_mass[0,0]) * unit_mass_in_cgs / Msun_in_cgs, linewidth=0.5, color='k', ms=0.5, label='swift')
#plot(time_Gyr,total_element_mass[:,0],linewidth=0.5,color='r',label='eagle test')
#xlabel("${\\rm{Time}} (Gyr)$", labelpad=0)
#ylabel("Change in element mass of gas particles (Msun)", labelpad=2)
#ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#legend()

savefig("box_evolution.png", dpi=200)




