###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
"""
Plot the output of testEOS to show how the equation of state pressure and sound
speed varies with density and specific internal energy.

Usage:
    python  testEOS.py  (mat_id)

    Sys args (optional):
        mat_id | int | Material ID, see equation_of_state.h for the options.
            Default: 201 (= HM80_ice).

Text file contains:
    header
    num_rho  num_u  mat_id                      # Header info
    rho_0   rho_1   rho_2   ...   rho_num_rho   # Array of densities, rho
    u_0     u_1     u_2     ...   u_num_u       # Array of energies, u
    P_0_0   P_0_1   ...     P_0_num_u           # Array of pressures, P(rho, u)
    P_1_0   ...     ...     P_1_num_u
    ...     ...     ...     ...
    P_num_rho_0     ...     P_num_rho_num_u
    c_0_0   c_0_1   ...     c_0_num_u           # Array of sound speeds, c(rho, u)
    c_1_0   ...     ...     c_1_num_u
    ...     ...     ...     ...
    c_num_rho_0     ...     c_num_rho_num_u

Note that the values tested extend beyond the range that most EOS are
designed for (e.g. outside table limits), to help test the EOS in cases of
unexpected particle behaviour.
"""

# ========
# Modules and constants
# ========

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

# Material types (copied from src/equation_of_state/planetary/equation_of_state.h)
type_factor = 100
Di_type = {"Til": 1, "HM80": 2, "SESAME": 3}
Di_material = {
    # Tillotson
    "Til_iron": Di_type["Til"] * type_factor,
    "Til_granite": Di_type["Til"] * type_factor + 1,
    "Til_water": Di_type["Til"] * type_factor + 2,
    # Hubbard & MacFarlane (1980) Uranus/Neptune
    "HM80_HHe": Di_type["HM80"] * type_factor,  # Hydrogen-helium atmosphere
    "HM80_ice": Di_type["HM80"] * type_factor + 1,  # H20-CH4-NH3 ice mix
    "HM80_rock": Di_type["HM80"] * type_factor + 2,  # SiO2-MgO-FeS-FeO rock mix
    # SESAME
    "SESAME_iron": Di_type["SESAME"] * type_factor,  # 2140
    "SESAME_basalt": Di_type["SESAME"] * type_factor + 1,  # 7530
    "SESAME_water": Di_type["SESAME"] * type_factor + 2,  # 7154
    "SS08_water": Di_type["SESAME"] * type_factor + 3,  # Senft & Stewart (2008)
}
# Invert so the mat_id are the keys
Di_mat_id = {mat_id: mat for mat, mat_id in list(Di_material.items())}

# Unit conversion
Ba_to_Mbar = 1e-12
erg_g_to_J_kg = 1e-4
cm_s_to_m_s = 1e-2

if __name__ == "__main__":
    # Sys args
    try:
        mat_id = int(sys.argv[1])
    except IndexError:
        mat_id = Di_material["HM80_ice"]

    # Check the material
    try:
        mat = Di_mat_id[mat_id]
        print(mat)
        sys.stdout.flush()
    except KeyError:
        print("Error: unknown material ID! mat_id = %d" % mat_id)
        print("Materials:")
        for mat_id, mat in sorted(Di_mat_id.items()):
            print("  %s%s%d" % (mat, (20 - len("%s" % mat)) * " ", mat_id))

    filename = "testEOS_rho_u_P_c_%d.txt" % mat_id

    # Load the header info and density and energy arrays
    with open(filename) as f:
        f.readline()
        num_rho, num_u, mat_id = np.array(f.readline().split(), dtype=int)
        A1_rho = np.array(f.readline().split(), dtype=float)
        A1_u = np.array(f.readline().split(), dtype=float)

    # Load pressure and soundspeed arrays
    A2_P = np.loadtxt(filename, skiprows=4)
    A2_c = A2_P[num_rho:]
    A2_P = A2_P[:num_rho]

    # Convert energies from cgs to SI
    A1_u *= erg_g_to_J_kg
    # Convert pressures from cgs (Barye) to Mbar
    A2_P *= Ba_to_Mbar
    # Convert sound speeds from cgs to SI
    A1_u *= cm_s_to_m_s

    # Check that the numbers are right
    assert A1_rho.shape == (num_rho,)
    assert A1_u.shape == (num_u,)
    assert A2_P.shape == (num_rho, num_u)
    assert A2_c.shape == (num_rho, num_u)

    # Plot
    # Pressure: P(rho) at fixed u
    plt.figure(figsize=(7, 7))
    ax = plt.gca()

    A1_colour = matplotlib.cm.rainbow(np.linspace(0, 1, num_u))

    for i_u, u in enumerate(A1_u):
        if i_u % 10 == 0:
            plt.plot(A1_rho, A2_P[:, i_u], c=A1_colour[i_u], label=r"%.2e" % u)
        else:
            plt.plot(A1_rho, A2_P[:, i_u], c=A1_colour[i_u])

    plt.legend(title="Sp. Int. Energy (J kg$^{-1}$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Density (g cm$^{-3}$)")
    plt.ylabel(r"Pressure (Mbar)")
    plt.title(mat)
    plt.tight_layout()

    plt.savefig("testEOS_P_%d.png" % mat_id)
    plt.close()

    # Sound speed: c(rho) at fixed u
    plt.figure(figsize=(7, 7))
    ax = plt.gca()

    A1_colour = matplotlib.cm.rainbow(np.linspace(0, 1, num_u))

    for i_u, u in enumerate(A1_u):
        if i_u % 10 == 0:
            plt.plot(A1_rho, A2_c[:, i_u], c=A1_colour[i_u], label=r"%.2e" % u)
        else:
            plt.plot(A1_rho, A2_c[:, i_u], c=A1_colour[i_u])

    plt.legend(title="Sp. Int. Energy (J kg$^{-1}$)")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Density (g cm$^{-3}$)")
    plt.ylabel(r"Sound Speed (m s^{-1})")
    plt.title(mat)
    plt.tight_layout()

    plt.savefig("testEOS_c_%d.png" % mat_id)
    plt.close()
