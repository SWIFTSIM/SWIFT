#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 15:45:33 2025

@author: roduit
"""

import swiftsimio as sw
import numpy as np
import matplotlib.pyplot as plt
import unyt
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid


filename = "snap/snapshot_0280.hdf5"
data = sw.load(filename)

m = data.gas.masses.to(unyt.Msun)
ids = data.gas.particle_ids
grad_Z = data.gas.metallicities_gradients
norm_grad_Z = np.linalg.norm(grad_Z, axis=1)

part = 18501


#%% 3D Scatter plot of the stars
# Create a 3D scatter plot
fig = plt.figure(figsize=(10, 8), num = 2)
ax = fig.add_subplot(111, projection='3d')

# Scatter plot
x = data.gas.coordinates[:, 0].value
y = data.gas.coordinates[:, 1].value
z = data.gas.coordinates[:, 2].value
a = data.gas.norm_diffusion_fluxes[:, 0].to(unyt.Msun/(unyt.Gyr*unyt.kpc**2)).value
d = np.log10(a)
# sc = ax.scatter(x, y, z, c=d, alpha=0.7, label='Stars')
sc = ax.scatter(x, y, z, c=a, alpha=0.7, label='Stars')

# Set labels and title
ax.set_xlabel('x [kpc]')
ax.set_ylabel('y [kpc]')
ax.set_zlabel('z [kpc]')
ax.set_title('3D Scatter Plot of Star Positions')

xmin = 0.092
xmax = 0.095
ymin = 0.429
ymax = 0.436
# xmin = 0.095
# xmax = 0.098
# ymin = 0.429
# ymax = 0.431

# ax.set_xlim([xmin, xmax])
# ax.set_ylim([ymin, ymax])

# Add colorbar
cbar = fig.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
cbar.set_label('Log Diffusion Flux [M$_\odot$/(Gyr kpc$^2)$]')


#%%
mask = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
x2 = data.gas.coordinates[:, 0].value[mask]
y2 = data.gas.coordinates[:, 1].value[mask]
z2 = data.gas.coordinates[:, 2].value[mask]
a2 = a[mask]
d2 = np.log10(a2)


fig = plt.figure(figsize=(10, 8), num = 3)
ax = fig.add_subplot(111, projection='3d')
# sc = ax.scatter(x2, y2, z2, c=d2, alpha=0.7, label='Stars')
sc = ax.scatter(x2, y2, z2, c=a2, alpha=0.7, label='Stars')


ax.set_xlabel('x [kpc]')
ax.set_ylabel('y [kpc]')
ax.set_zlabel('z [kpc]')

# Add colorbar
cbar = fig.colorbar(sc, ax=ax, shrink=0.5, aspect=10)
cbar.set_label('Log Diffusion Flux [M$_\odot$/(Gyr kpc$^2)$]')

#%%

import re

file = "/Users/roduit/swiftsim/examples/ChemistryTests/HomogeneousBox/output.log"
# file = "/Users/roduit/swiftsim/examples/ChemistryTests/HomogeneousBox/data.log"


# Dictionary to store data grouped by the first number
grouped_data = {}

# Open the file for reading
with open(file, 'r') as file:
    for line in file:
        # Check if the line matches the pattern of the unwanted line (numerical data without keywords)
        if re.match(r'^\s*\d+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+(\d+\s+){6}[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+$', line):
            continue  # Skip this line as it matches the unwanted format
        
        # Check if the line contains the relevant keywords (e.g., 'mindt', 'totflux', 'mass_flux', etc.)
        if "mindt" in line and "totflux" in line and "mass_flux" in line:
            # Use a regular expression to capture the relevant data
            match = re.search(r'\[\d+\.\d+\]\s+runner_iact_chemistry_fluxes_common: \[(\d+)\s+\|\s+(\d+),\s+(\d+)\]\s+mindt\s*=\s*([0-9.e+-]+),\s*totflux\s*=\s*([0-9.e+-]+),\s*mass_flux\s*=\s*([0-9.e+-]+),\s*mode\s*=\s*(\d+),\s*chi->flux_dt\s*=\s*([0-9.e+-]+),\s*chj->flux_dt\s*=\s*([0-9.e+-]+),\s*m_Zi\s*=\s*([0-9.e+-]+),\s*m_Zj\s*=\s*([0-9.e+-]+),\s*mi\s*=\s*([0-9.e+-]+),\s*mj\s*=\s*([0-9.e+-]+)', line)
            
            if match:
                # Extract the time, id_i, id_j and other parameters
                time = int(match.group(1))  # Extract time (first number in the second brackets)
                id_i = int(match.group(2))  # Extract id_i (second number in the second brackets)
                id_j = int(match.group(3))  # Extract id_j (third number in the second brackets)
                mindt = float(match.group(4))
                totflux = float(match.group(5))
                mass_flux = float(match.group(6))
                mode = int(match.group(7))
                chi_flux_dt = float(match.group(8))
                chj_flux_dt = float(match.group(9))
                m_Zi = float(match.group(10))
                m_Zj = float(match.group(11))
                mi = float(match.group(12))
                mj = float(match.group(13))
                
                # Group data by the extracted time
                if time not in grouped_data:
                    grouped_data[time] = []

                # Append the extracted data to the corresponding group
                grouped_data[time].append({
                    'id_i': id_i,
                    'id_j': id_j,
                    'mindt': mindt,
                    'totflux': totflux,
                    'mass_flux': mass_flux,
                    'mode': mode,
                    'chi_flux_dt': chi_flux_dt,
                    'chj_flux_dt': chj_flux_dt,
                    'm_Zi': m_Zi,
                    'm_Zj': m_Zj,
                    'mi': mi,
                    'mj': mj
                })

#%%

file = "/Users/roduit/swiftsim/examples/ChemistryTests/HomogeneousBox/output.log"

# Dictionary to store data grouped by the first number
grouped_by_id_i = {}
grouped_by_id_j = {}

# Open the file for reading
with open(file, 'r') as file:
    for line in file:
        # Check if the line matches the pattern of the unwanted line (numerical data without keywords)
        if re.match(r'^\s*\d+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+(\d+\s+){6}[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+$', line):
            continue  # Skip this line as it matches the unwanted format
        
        # Check if the line contains the relevant keywords (e.g., 'mindt', 'totflux', 'mass_flux', etc.)
        if "mindt" in line and "totflux" in line and "mass_flux" in line:
            # Use a regular expression to capture the relevant data
            match = re.search(r'\[\d+\.\d+\]\s+runner_iact_chemistry_fluxes_common: \[(\d+)\s+\|\s+(\d+),\s+(\d+)\]\s+mindt\s*=\s*([0-9.e+-]+),\s*totflux\s*=\s*([0-9.e+-]+),\s*mass_flux\s*=\s*([0-9.e+-]+),\s*mode\s*=\s*(\d+),\s*chi->flux_dt\s*=\s*([0-9.e+-]+),\s*chj->flux_dt\s*=\s*([0-9.e+-]+),\s*m_Zi\s*=\s*([0-9.e+-]+),\s*m_Zj\s*=\s*([0-9.e+-]+),\s*mi\s*=\s*([0-9.e+-]+),\s*mj\s*=\s*([0-9.e+-]+)', line)
            
            if match:
                # Extract the time, id_i, id_j, and other parameters
                time = int(match.group(1))  # Extract time (first number in the second brackets)
                id_i = int(match.group(2))  # Extract id_i (second number in the second brackets)
                id_j = int(match.group(3))  # Extract id_j (third number in the second brackets)
                mindt = float(match.group(4))
                totflux = float(match.group(5))
                mass_flux = float(match.group(6))
                mode = int(match.group(7))
                chi_flux_dt = float(match.group(8))
                chj_flux_dt = float(match.group(9))
                m_Zi = float(match.group(10))
                m_Zj = float(match.group(11))
                mi = float(match.group(12))
                mj = float(match.group(13))
                
                # Group data by id_i (particle i)
                if id_i not in grouped_by_id_i:
                    grouped_by_id_i[id_i] = []
                grouped_by_id_i[id_i].append({
                    'time': time,
                    'id_j': id_j,
                    'mindt': mindt,
                    'totflux': totflux,
                    'mass_flux': mass_flux,
                    'mode': mode,
                    'chi_flux_dt': chi_flux_dt,
                    'chj_flux_dt': chj_flux_dt,
                    'm_Zi': m_Zi,
                    'm_Zj': m_Zj,
                    'mi': mi,
                    'mj': mj
                })
                
                # Group data by id_j (particle j)
                if id_j not in grouped_by_id_j:
                    grouped_by_id_j[id_j] = []
                grouped_by_id_j[id_j].append({
                    'time': time,
                    'id_i': id_i,
                    'mindt': mindt,
                    'totflux': totflux,
                    'mass_flux': mass_flux,
                    'mode': mode,
                    'chi_flux_dt': chi_flux_dt,
                    'chj_flux_dt': chj_flux_dt,
                    'm_Zi': m_Zi,
                    'm_Zj': m_Zj,
                    'mi': mi,
                    'mj': mj
                })


