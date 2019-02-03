###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2019 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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

# Plot the snapshots from the example giant impact on Uranus, showing the 
# particles in a thin slice near z=0, coloured by their material, similarish
# (but not identical) to Fig. 2 in Kegerreis et al. (2018).

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
# import swiftsimio as sw
import h5py

font_size   = 20
params      = { 
    'axes.labelsize'    : font_size,
    'font.size'         : font_size,
    'xtick.labelsize'   : font_size,
    'ytick.labelsize'   : font_size,
    'text.usetex'       : True,
    'font.family'       : 'serif',
    }
matplotlib.rcParams.update(params)

# Snapshot output times
output_list = [4000, 9000, 14000, 20000, 30000, 40000]

# Material IDs ( = type_id * type_factor + unit_id )
type_factor = 100
type_HM80   = 2
id_body     = 10000
# Name and ID
Di_mat_id   = {
    'HM80_HHe'      : type_HM80 * type_factor,      # Hydrogen-helium atmosphere
    'HM80_ice'      : type_HM80 * type_factor + 1,  # H20-CH4-NH3 ice mix
    'HM80_ice_2'    : type_HM80 * type_factor + 1 + id_body,
    'HM80_rock'     : type_HM80 * type_factor + 2,  # SiO2-MgO-FeS-FeO rock mix
    'HM80_rock_2'   : type_HM80 * type_factor + 2 + id_body,
    }
# ID and colour
Di_id_colour    = {
    Di_mat_id['HM80_HHe']       : '#33DDFF',
    Di_mat_id['HM80_ice']       : 'lightsteelblue',
    Di_mat_id['HM80_ice_2']     : '#A080D0',
    Di_mat_id['HM80_rock']      : 'slategrey',
    Di_mat_id['HM80_rock_2']    : '#706050',
    }
   
def get_snapshot_slice(snapshot):
    """ Load and select the particles to plot. """
    # Load particle data
    # data    = load("uranus_1e6_%06d.hdf5" % snapshot)
    # id      = data.gas.particle_ids
    # pos     = data.gas.coordinates
    # mat_id  = data.gas.material
    with h5py.File("uranus_1e6_%06d.hdf5" % snapshot, 'r') as f:
        id      = f['PartType0/ParticleIDs'].value
        pos     = (f['PartType0/Coordinates'].value
                   - 0.5 * f['Header'].attrs['BoxSize'])
        mat_id  = f['PartType0/MaterialID'].value
                   
    # Edit the material ID of particles in the impactor
    num_in_target   = 869104
    sel_id          = np.where(num_in_target < id)[0]
    mat_id[sel_id]  += id_body
    
    # Select particles in a thin slice around z=0
    z_min   = -0.1
    z_max   = 0.1
    sel_z   = np.where((z_min < pos[:, 2]) & (pos[:, 2] < z_max))[0]
    pos     = pos[sel_z]
    mat_id  = mat_id[sel_z]
    
    return pos, mat_id

def plot_snapshot_slice(pos, mat_id):
    """ Plot the particles, coloured by their material. """
    colour  = np.empty(len(pos), dtype=object)
    for id, c in Di_id_colour.items():
        sel_c           = np.where(mat_id == id)[0]
        colour[sel_c]   = c

    ax.scatter(pos[:, 0], pos[:, 1], c=colour, edgecolors='none', marker='.', 
               s=10, alpha=0.5, zorder=0)

# Set up the figure
fig     = plt.figure(figsize=(12, 8))
gs      = matplotlib.gridspec.GridSpec(2, 3)
axes    = [plt.subplot(gs[i_y, i_x]) for i_y in range(2) for i_x in range(3)]

# Plot each snapshot
for i_ax, ax in enumerate(axes):
    plt.sca(ax)
    ax.set_rasterization_zorder(1)
    
    # Load and select the particles to plot
    pos, mat_id = get_snapshot_slice(output_list[i_ax])
    
    # Plot the particles, coloured by their material
    plot_snapshot_slice(pos, mat_id)
    
    # Axes etc.
    ax.set_aspect('equal')
    ax.set_facecolor('k')
    
    ax.set_xlim(-13, 13)
    ax.set_ylim(-13, 13)

    if i_ax in [0, 3]:
        ax.set_ylabel(r"y Postion $(R_\oplus)$")
    else: 
        ax.set_yticklabels([])
    if 2 < i_ax:
        ax.set_xlabel(r"x Postion $(R_\oplus)$")
    else: 
        ax.set_xticklabels([])
        
    # Corner time labels 
    x   = ax.get_xlim()[0] + 0.04 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    y   = ax.get_ylim()[0] + 0.89 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.text(x, y, "%.1f h" % (output_list[i_ax] / 60**2), color='w')

plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout()

# Save
plt.savefig("uranus_1e6.pdf", dpi=200)