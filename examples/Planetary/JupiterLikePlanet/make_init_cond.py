################################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
# 		        2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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

"""Create initial conditions for settling, using WoMa. See README.md for more info."""

import numpy as np
import h5py
import woma

# Number of particles
N = 10 ** 7
N_label = "n%d" % (10 * np.log10(N))

# Earth units
M_E = 5.9724e24  # kg
R_E = 6.3710e6  # m

# Set profile inputs
M_t = 308 * M_E
target_prof = woma.Planet(
    name="target",
    A1_mat_layer=["AQUA", "CD21_HHe"],
    A1_T_rho_type=["adiabatic", "adiabatic"],
    M=M_t,
    A1_M_layer=[10 * M_E, 298 * M_E],
    P_s=1e5,
    T_s=165,
    num_prof=10000,
)

# Load material tables
woma.load_eos_tables(np.unique(target_prof.A1_mat_layer))

# Compute profiles
target_prof.gen_prof_L2_find_R_R1_given_M1_M2(R_min=10.7 * R_E, R_max=11.2 * R_E)

# Save profile data
target_prof.save("demo_target_profile.hdf5")

# Place particles
target = woma.ParticlePlanet(target_prof, N, seed=12345)

print()
print("N_target     = %d" % target.N_particles)

# Save the settling initial conditions
file_to_SI = woma.Conversions(m=1e24, l=1e6, t=1)
target.save(
    "demo_target_%s.hdf5" % N_label,
    boxsize=30 * R_E,
    file_to_SI=file_to_SI,
    do_entropies=True,
)
