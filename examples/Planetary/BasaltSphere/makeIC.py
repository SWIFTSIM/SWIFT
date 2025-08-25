###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

import h5py
from numpy import *
import random

# Generates a swift IC file
fileOutputName = "impact.hdf5"
# ---------------------------------------------------

# from AKIKO NAKAMURA AND AKIRA FUJIWARA 1991
L_target = 300#200  # Number of particles on side of grid
boxsize = 20 / 100
vol = boxsize**3
numPart_grid = L_target * L_target * L_target

target_radius = 3 / 100 #cm
rho_target = 2700
mat_id_target = 1000
u_target = 0.07407407408291543#0.07407407408291543#3e5

# Place particles in grid
L_grid = L_target
A2_pos_grid = zeros((numPart_grid, 3))
for i in range(L_grid):
    for j in range(L_grid):
         for k in range(L_grid):
            index = i * L_grid * L_grid + j * L_grid + k
            A2_pos_grid[index, 0] = (i / (float(L_grid))) * boxsize
            A2_pos_grid[index, 1] = (j / (float(L_grid))) * boxsize
            A2_pos_grid[index, 2] = (k / (float(L_grid))) * boxsize


# Select sphere 
A1_select_target = (A2_pos_grid[:, 0] - boxsize/2)**2 + (A2_pos_grid[:, 1] - boxsize/2)**2 + (A2_pos_grid[:, 2] - boxsize/2)**2 <= target_radius**2


A2_pos_target = A2_pos_grid[A1_select_target, :]

num_target = len(A2_pos_target)
A2_vel_target = zeros((num_target, 3))
A1_h_target = ones(num_target) * (boxsize / L_grid) * 1.2348
A1_m_target = ones(num_target) * vol * rho_target / numPart_grid
A1_mat_id_target = ones(num_target) * mat_id_target
A1_rho_target = ones(num_target) * rho_target
A1_u_target = ones(num_target) * u_target

print(num_target)


# Flaw parameters from B&A 1994
m_flaw_target = 8#8.5
k_flaw_target = 5e28 * 100**3 # cm^-3 => m^-3
Volume = num_target * vol / numPart_grid
# a randomly chosen fixed seed for testing
seed = 1#12
random.seed(seed)
hardcoded_max_flaws = 40 # log(1e9)


all_filled = 0
N_flaw = 1

A1_particle_flaws_target = zeros(num_target)
A2_activation_thresholds_target = zeros((num_target, hardcoded_max_flaws))

# Generate activation_thresholds until every particle has a flaw
while not all_filled:
    # Generates random index in sequence set by seed
    idx = random.randrange(0, num_target)
    
    current_activation_threshold = (N_flaw / (k_flaw_target * Volume))**(1 / m_flaw_target)
    
    A2_activation_thresholds_target[idx, int(A1_particle_flaws_target[idx])] = current_activation_threshold
    
    A1_particle_flaws_target[idx] +=1 
    N_flaw += 1
    if 0 not in A1_particle_flaws_target:
        all_filled = 1
    
# Print out number of flaws in particle with the most flaws as a sanity check   
max_flaws = int(max(A1_particle_flaws_target))
print(max_flaws)
    

# Now make impactor in same way as target

rho_impactor = 1180#1800
impactor_mass = 0.2 / 1000 #g => kg
impactor_radius = cbrt(impactor_mass / (rho_impactor) / (4/3 * pi))#0.35 / 100 #cm 

#rho_impactor = rho_target
mat_id_impactor = 1001
u_impactor = 0.07407407408291543#1655347#0.07407407408291543#3e5

L_impactor = L_target#int(L_target * cbrt(rho_impactor / rho_target))
numPart_grid = L_impactor * L_impactor * L_impactor
L_grid = L_impactor
A2_pos_grid = zeros((numPart_grid, 3))
for i in range(L_grid):
    for j in range(L_grid):
         for k in range(L_grid):
            index = i * L_grid * L_grid + j * L_grid + k
            A2_pos_grid[index, 0] = (i / (float(L_grid))) * boxsize
            A2_pos_grid[index, 1] = (j / (float(L_grid))) * boxsize
            A2_pos_grid[index, 2] = (k / (float(L_grid))) * boxsize

A1_select_impactor = (A2_pos_grid[:, 0] - boxsize/2)**2 + (A2_pos_grid[:, 1] - boxsize/2)**2 + (A2_pos_grid[:, 2] - boxsize/2)**2 <= impactor_radius**2

A2_pos_impactor = A2_pos_grid[A1_select_impactor, :]

num_impactor = len(A2_pos_impactor)
A2_vel_impactor = zeros((num_impactor, 3))
A1_h_impactor = ones(num_impactor) * (boxsize / L_grid) * 1.2348
A1_m_impactor = ones(num_impactor) * vol * rho_impactor / numPart_grid
A1_mat_id_impactor = ones(num_impactor) * mat_id_impactor
A1_rho_impactor = ones(num_impactor) * rho_impactor
A1_u_impactor = ones(num_impactor) * u_impactor

print(num_impactor)

vx = 3.2e5 / 100#3.2e5 / 100 #cm/s
A2_vel_impactor[:,0] = vx


# Set initial offset in x direction of impactor fom target 
initial_offset = 0.2 * boxsize
A2_pos_impactor[:,0] -= initial_offset

# Set impact angle
angle = 30
impact_parameter_b = (target_radius + impactor_radius) * sin(deg2rad(angle))
A2_pos_impactor[:, 1] += impact_parameter_b

# Impactor has no flaws
A1_particle_flaws_impactor = zeros(num_impactor)
A2_activation_thresholds_impactor = zeros((num_impactor, hardcoded_max_flaws))

pos = append(A2_pos_target, A2_pos_impactor, axis=0)
vel = append(A2_vel_target, A2_vel_impactor, axis=0)
h = append(A1_h_target, A1_h_impactor, axis=0)
m = append(A1_m_target, A1_m_impactor)
mat = append(A1_mat_id_target, A1_mat_id_impactor)
rho = append(A1_rho_target, A1_rho_impactor)
u = append(A1_u_target, A1_u_impactor)

particle_flaws = append(A1_particle_flaws_target, A1_particle_flaws_impactor)
activation_thresholds = append(A2_activation_thresholds_target, A2_activation_thresholds_impactor, axis=0)

numPart = len(h)
ids = linspace(1, numPart, numPart)

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [boxsize, boxsize, boxsize]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 100.
grp.attrs["Unit mass in cgs (U_M)"] = 1000.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = fileOutput.create_group("/PartType0")
ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
ds[()] = pos
ds = grp.create_dataset("Velocities", (numPart, 3), "f")
ds[()] = vel
ds = grp.create_dataset("Masses", (numPart, 1), "f")
ds[()] = m.reshape((numPart, 1))
ds = grp.create_dataset("Density", (numPart, 1), "f")
ds[()] = rho.reshape((numPart, 1))
ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
ds[()] = h.reshape((numPart, 1))
ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
ds[()] = u.reshape((numPart, 1))
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids.reshape((numPart, 1))
ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
ds[()] = mat.reshape((numPart, 1))

# Save number of flaws
ds = grp.create_dataset("NumFlaws", (numPart, 1), "i")
ds[()] = particle_flaws.reshape((numPart, 1))
# Save activation thresholds
ds = grp.create_dataset("ActivationThresholds", (numPart, hardcoded_max_flaws), "f")
ds[()] = activation_thresholds


fileOutput.close()
