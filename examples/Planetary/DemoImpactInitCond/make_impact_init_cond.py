################################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
################################################################################

"""Create the initial conditions for the DemoImpact example. See README.md for more info."""

import numpy as np
import h5py
import woma

# Number of particles
N = 10**6
N_label = "n%d" % (10 * np.log10(N))

# Earth units
M_E = 5.9724e24  # kg
R_E = 6.3710e6  # m

# Set profile inputs
M_t = 0.887 * M_E
M_i = 0.133 * M_E

# Load the settled particle planets
def load_snapshot_data(filename):
    with h5py.File(filename, "r") as f:
        # Units from file metadata
        file_to_SI = woma.Conversions(
            m=float(f["Units"].attrs["Unit mass in cgs (U_M)"]) * 1e-3,
            l=float(f["Units"].attrs["Unit length in cgs (U_L)"]) * 1e-2,
            t=float(f["Units"].attrs["Unit time in cgs (U_t)"]),
        )

        # Particle data (converted to SI)
        A2_pos = (
            np.array(f["PartType0/Coordinates"][()])
            - 0.5 * f["Header"].attrs["BoxSize"]
        ) * file_to_SI.l
        A1_m = np.array(f["PartType0/Masses"][()]) * file_to_SI.m
        A1_h = np.array(f["PartType0/SmoothingLengths"][()]) * file_to_SI.l
        A1_rho = np.array(f["PartType0/Densities"][()]) * file_to_SI.rho
        A1_P = np.array(f["PartType0/Pressures"][()]) * file_to_SI.P
        A1_u = np.array(f["PartType0/InternalEnergies"][()]) * file_to_SI.u
        A1_mat_id = np.array(f["PartType0/MaterialIDs"][()])

        return A2_pos, A1_m, A1_h, A1_rho, A1_P, A1_u, A1_mat_id


snapshot_id = 5
A2_pos_t, A1_m_t, A1_h_t, A1_rho_t, A1_P_t, A1_u_t, A1_mat_id_t = load_snapshot_data(
    "snapshots/demo_target_%s_%04d.hdf5" % (N_label, snapshot_id)
)
A2_pos_i, A1_m_i, A1_h_i, A1_rho_i, A1_P_i, A1_u_i, A1_mat_id_i = load_snapshot_data(
    "snapshots/demo_impactor_%s_%04d.hdf5" % (N_label, snapshot_id)
)

# Impact initial conditions (target rest frame)
# Collide at 45 degrees, at the mutual escape speed, start 1 hour before contact
A1_pos_t, A1_vel_t = np.zeros(3), np.zeros(3)
A1_pos_i, A1_vel_i = woma.impact_pos_vel_b_v_c_t(
    b=np.sin(45 * np.pi / 180),
    v_c=1,
    units_v_c="v_esc",
    t=3600,
    R_t=1.000 * R_E,
    R_i=0.566 * R_E,
    M_t=0.887 * M_E,
    M_i=0.133 * M_E,
)

# Shift to centre-of-mass frame
A1_pos_com = (M_t * A1_pos_t + M_i * A1_pos_i) / (M_t + M_i)
A1_vel_com = (M_t * A1_vel_t + M_i * A1_vel_i) / (M_t + M_i)

A1_pos_t -= A1_pos_com
A1_vel_t -= A1_vel_com
A1_pos_i -= A1_pos_com
A1_vel_i -= A1_vel_com

# Update particle positions
A2_pos_t[:] += A1_pos_t
A2_vel_t = np.zeros_like(A2_pos_t) + A1_vel_t
A2_pos_i[:] += A1_pos_i
A2_vel_i = np.zeros_like(A2_pos_i) + A1_vel_i

# Combine and save the particle data (label by the resolution)
with h5py.File("demo_impact_%s.hdf5" % N_label, "w") as f:
    woma.save_particle_data(
        f,
        A2_pos=np.append(A2_pos_t, A2_pos_i, axis=0),
        A2_vel=np.append(A2_vel_t, A2_vel_i, axis=0),
        A1_m=np.append(A1_m_t, A1_m_i),
        A1_h=np.append(A1_h_t, A1_h_i),
        A1_rho=np.append(A1_rho_t, A1_rho_i),
        A1_P=np.append(A1_P_t, A1_P_i),
        A1_u=np.append(A1_u_t, A1_u_i),
        A1_mat_id=np.append(A1_mat_id_t, A1_mat_id_i),
        boxsize=100 * R_E,
        file_to_SI=woma.Conversions(m=1e24, l=1e6, t=1),
    )
