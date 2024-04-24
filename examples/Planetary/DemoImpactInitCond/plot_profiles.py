###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
##############################################################################

"""Plot the post-settling planetary profiles of the DemoImpactInitCond simulations.

Note that, for standard SPH hydro schemes, especially at low resolution, the
standard issues that arise at discontinuities in material and density lead the
SPH particles to settle at slightly different densities near discontinuties.
For more info and to explore ways to resolve these issues, check out e.g.
Sandnes et al. (2024), Ruiz-Bonilla el al. (2022), and Kegerreis et al. (2019).
The overall profile of the settled SPH planet should still align with the input.
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import woma

# Number of particles
N = 10**6
N_label = "n%d" % (10 * np.log10(N))

# Earth units
M_E = 5.9724e24  # kg
R_E = 6.3710e6  # m

# Plotting options
font_size = 20
params = {
    "axes.labelsize": font_size,
    "font.size": font_size,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    "font.family": "serif",
}
matplotlib.rcParams.update(params)


def plot_profile_and_particles(profile, A1_r, A1_rho):
    """Plot the particles."""
    plt.figure(figsize=(7, 7))
    ax = plt.gca()

    # Profile
    ax.plot(profile.A1_r / R_E, profile.A1_rho)

    # Particles
    ax.scatter(A1_r / R_E, A1_rho, c="k", marker=".", s=1**2)

    ax.set_xlim(0, None)
    ax.set_xlabel(r"Radial distance ($R_\oplus$)")
    ax.set_ylabel(r"Density (kg m$^{-3}$)")

    plt.tight_layout()


if __name__ == "__main__":
    # Plot each snapshot
    for body in ["target", "impactor"]:
        # Load profiles
        profile = woma.Planet(load_file="demo_%s_profile.hdf5" % body)

        # Load the data
        snapshot_id = 5
        filename = "snapshots/demo_%s_%s_%04d.hdf5" % (body, N_label, snapshot_id)
        with h5py.File(filename, "r") as f:
            # Units from file metadata
            file_to_SI = woma.Conversions(
                m=float(f["Units"].attrs["Unit mass in cgs (U_M)"]) * 1e-3,
                l=float(f["Units"].attrs["Unit length in cgs (U_L)"]) * 1e-2,
                t=float(f["Units"].attrs["Unit time in cgs (U_t)"]),
            )

            # Particle data
            A2_pos = (
                np.array(f["PartType0/Coordinates"][()])
                - 0.5 * f["Header"].attrs["BoxSize"]
            ) * file_to_SI.l
            A1_r = np.sqrt(np.sum(A2_pos**2, axis=1))
            A1_rho = np.array(f["PartType0/Densities"][()]) * file_to_SI.rho

        # Plot the data
        plot_profile_and_particles(profile, A1_r, A1_rho)

        # Save the figure
        save = "demo_%s_%s_%04d_prof.png" % (body, N_label, snapshot_id)
        plt.savefig(save, dpi=200)
        plt.close()

        print("\rSaved %s" % save)
