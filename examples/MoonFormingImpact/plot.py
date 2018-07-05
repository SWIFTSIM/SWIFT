"""
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
###############################################################################

Plotting script for the Canonical Moon-Forming Giant Impact example.

Save a figure for each snapshot in `./plots/` then make them into a simple
animation with ffmpeg in `./`.

Usage:
    `$ python  plot.py  time_end  delta_time`

    Sys args:
        + `time_end` | (opt) int | The time of the last snapshot to plot.
            Default = 100000
        + `delta_time` | (opt) int | The time between successive snapshots.
            Default = 100
"""

from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import sys
import subprocess

# Particle array fields
dtype_picle = [
    ('m', float), ('x', float), ('y', float), ('z', float), ('v_x', float),
    ('v_y', float), ('v_z', float), ('ID', int), ('rho', float), ('u', float),
    ('phi', float), ('P', float), ('h', float), ('mat_ID', int), ('r', float)
    ]

s_to_hour   = 1 / 60**2

# Snapshot info
file_snap   = "./snapshots/moon_forming_impact_"
file_plot   = "./plots/moon_forming_impact_"
# Number of particles in the target body
num_target  = 9496

# Material types (copied from src/equation_of_state/planetary/equation_of_state.h)
type_factor = 100
Di_type = {
    'Til'       : 1,
    'HM80'      : 2,
    'ANEOS'     : 3,
    'SESAME'    : 4,
}
Di_material = {
    # Tillotson
    'Til_iron'      : Di_type['Til']*type_factor,
    'Til_granite'   : Di_type['Til']*type_factor + 1,
    'Til_water'     : Di_type['Til']*type_factor + 2,
    # Hubbard & MacFarlane (1980) Uranus/Neptune
    'HM80_HHe'      : Di_type['HM80']*type_factor,      # Hydrogen-helium atmosphere
    'HM80_ice'      : Di_type['HM80']*type_factor + 1,  # H20-CH4-NH3 ice mix
    'HM80_rock'     : Di_type['HM80']*type_factor + 2,  # SiO2-MgO-FeS-FeO rock mix
    # ANEOS
    'ANEOS_iron'        : Di_type['ANEOS']*type_factor,
    'MANEOS_forsterite' : Di_type['ANEOS']*type_factor + 1,
    # SESAME
    'SESAME_iron'   : Di_type['SESAME']*type_factor,
}

# Material offset for impactor particles
ID_imp  = 10000
# Material colours
Di_mat_colour = {
    # Target
    Di_material['Til_iron']             : 'orange',
    Di_material['Til_granite']          : '#FFF0E0',
    # Impactor
    Di_material['Til_iron'] + ID_imp    : 'dodgerblue',
    Di_material['Til_granite'] + ID_imp : '#A080D0',
    }


def load_snapshot(filename):
    """ Load the hdf5 snapshot file and return the structured particle array.
    """
    # Add extension if needed
    if (filename[-5:] != ".hdf5"):
        filename += ".hdf5"

    # Load the hdf5 file
    with h5py.File(filename, 'r') as f:
        header      = f['Header'].attrs
        A2_pos      = f['PartType0/Coordinates'].value
        A2_vel      = f['PartType0/Velocities'].value

        # Structured array of all particle data
        A2_picle    = np.empty(header['NumPart_Total'][0],
                               dtype=dtype_picle)

        A2_picle['x']       = A2_pos[:, 0]
        A2_picle['y']       = A2_pos[:, 1]
        A2_picle['z']       = A2_pos[:, 2]
        A2_picle['v_x']     = A2_vel[:, 0]
        A2_picle['v_y']     = A2_vel[:, 1]
        A2_picle['v_z']     = A2_vel[:, 2]
        A2_picle['m']       = f['PartType0/Masses'].value
        A2_picle['ID']      = f['PartType0/ParticleIDs'].value
        A2_picle['rho']     = f['PartType0/Density'].value
        A2_picle['u']       = f['PartType0/InternalEnergy'].value
        A2_picle['phi']     = f['PartType0/Potential'].value
        A2_picle['P']       = f['PartType0/Pressure'].value
        A2_picle['h']       = f['PartType0/SmoothingLength'].value
        A2_picle['mat_ID']  = f['PartType0/MaterialID'].value

    return A2_picle


def process_particles(A2_picle, num_target):
    """ Modify things like particle units, material IDs, and coordinate origins.
    """
    # Offset material IDs for impactor particles
    A2_picle['mat_ID'][A2_picle['ID'] >= num_target] += ID_imp

    # Shift coordinates to the centre of the target's core's mass and momentum
    sel_tar  = np.where(A2_picle['mat_ID'] == Di_material['Til_iron'])[0]

    # Centre of mass
    m_tot   = np.sum(A2_picle[sel_tar]['m'])
    x_com   = np.sum(A2_picle[sel_tar]['m'] * A2_picle[sel_tar]['x']) / m_tot
    y_com   = np.sum(A2_picle[sel_tar]['m'] * A2_picle[sel_tar]['y']) / m_tot
    z_com   = np.sum(A2_picle[sel_tar]['m'] * A2_picle[sel_tar]['z']) / m_tot

    # Change origin to the centre-of-mass
    A2_picle['x']   -= x_com
    A2_picle['y']   -= y_com
    A2_picle['z']   -= z_com
    A2_picle['r']   = np.sqrt(
        A2_picle['x']**2 + A2_picle['y']**2 + A2_picle['z']**2
        )

    # Centre of momentum
    v_x_com = np.sum(A2_picle[sel_tar]['m'] * A2_picle[sel_tar]['v_x']) / m_tot
    v_y_com = np.sum(A2_picle[sel_tar]['m'] * A2_picle[sel_tar]['v_y']) / m_tot
    v_z_com = np.sum(A2_picle[sel_tar]['m'] * A2_picle[sel_tar]['v_z']) / m_tot

    # Change to the centre-of-momentum frame of reference
    A2_picle['v_x'] -= v_x_com
    A2_picle['v_y'] -= v_y_com
    A2_picle['v_z'] -= v_z_com

    return A2_picle


def plot_snapshot(A2_picle, filename, time, ax_lim=100, dz=0.1):
    """ Plot the snapshot particles and save the figure.
    """
    # Add extension if needed
    if (filename[-5:] != ".png"):
        filename += ".png"

    fig = plt.figure(figsize=(9, 9))
    ax  = fig.add_subplot(111, aspect='equal')

    # Plot slices in z below zero
    for z in np.arange(-ax_lim, 0, dz):
        sel_z       = np.where((z < A2_picle['z']) & (A2_picle['z'] < z+dz))[0]
        A2_picle_z  = A2_picle[sel_z]

        # Plot each material
        for mat_ID, colour in Di_mat_colour.iteritems():
            sel_col = np.where(A2_picle_z['mat_ID'] == mat_ID)[0]

            ax.scatter(
                A2_picle_z[sel_col]['x'], A2_picle_z[sel_col]['y'],
                c=colour, edgecolors='none', marker='.', s=50, alpha=0.7
                )

    # Axes etc.
    ax.set_axis_bgcolor('k')

    ax.set_xlabel("x Position ($R_\oplus$)")
    ax.set_ylabel("y Position ($R_\oplus$)")

    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)

    plt.text(
        -0.92*ax_lim, 0.85*ax_lim, "%.1f h" % (time*s_to_hour), fontsize=20,
        color='w'
        )

    # Font sizes
    for item in (
        [ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() +
        ax.get_yticklabels()
        ):
        item.set_fontsize(20)

    plt.tight_layout()

    plt.savefig(filename)
    plt.close()


if __name__ == '__main__':
    # Sys args
    try:
        time_end    = int(sys.argv[1])

        try:
            delta_time  = int(sys.argv[2])
        except IndexError:
            delta_time  = 100
    except IndexError:
        time_end    = 100000
        delta_time  = 100

    # Load and plot each snapshot
    for i_snap in range(int(time_end/delta_time) + 1):
        snap_time   = i_snap * delta_time
        print "\rPlotting snapshot %06d (%d of %d)" % (
            snap_time, i_snap+1, int(time_end/delta_time)
            ),
        sys.stdout.flush()

        # Load particle data
        filename    = "%s%06d" % (file_snap, snap_time)
        A2_picle    = load_snapshot(filename)

        # Process particle data
        A2_picle    = process_particles(A2_picle, num_target)

        # Plot particles
        filename    = "%s%06d" % (file_plot, snap_time)
        plot_snapshot(A2_picle, filename, snap_time)

    # Animation
    command = (
        "ffmpeg -framerate 12 -i plots/moon_forming_impact_%*.png -r 25 "
        "anim.mpg -y"
        )
    print "\n%s\n" % command
    subprocess.check_output(command, shell=True)





























