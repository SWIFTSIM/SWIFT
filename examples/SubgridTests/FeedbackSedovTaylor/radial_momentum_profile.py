###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
import argparse
import os
import swiftsimio as sw
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import unyt as u

#%%

def list_snapshots(folder_path):
    """Lists all snapshot files in the folder."""
    snapshots = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.startswith("snapshot_") and f.endswith(".hdf5")
    ]
    return sorted(snapshots)


#%%

class ExplicitDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Formatter that prints default arguments if provided."""

    def _get_help_string(self, action):
        if isinstance(action.default, np.ndarray) or action.default not in (None, False):
            return super()._get_help_string(action)
        return action.help


class RawTextArgumentDefaultsHelpFormatter(ExplicitDefaultsHelpFormatter,
                                           argparse.RawTextHelpFormatter):
    """Combines raw text formatting with default argument printing."""
    pass



def parse_option():
    description = """"
Compute the metallicity distribution function (MDF).
"""
    epilog = """
Examples:
--------
madr_plot_mdf.py UFD_0172.hdf5 --x_min -20 --x_max 0 --output_location "MDF" --n_bins 100

"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=RawTextArgumentDefaultsHelpFormatter)

    parser.add_argument("files",
                        nargs="+",
                        type=str,
                        help="File name(s).")

    parser.add_argument("-o",
                        action="store",
                        type=str,
                        dest="output_filename",
                        default=None,
                        help="Name of the output file. Use it if you only give ONE input file. Otherwise, "
                        "all files will have the same name.")

    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File not found: {f}")

    return args, files


#%%
# Parse arguments
args, files = parse_option()

# Use output location and filename arguments, with defaults
output_filename = args.output_filename if args.output_filename else "momentum_profile.png"

# Parent directory containing multiple runs
base_folders = ["./"]

all_simulations_momentum = []

for base_folder in base_folders:
    # For each file in the provided list of files
    for file in files:
        snapshot = file
        momentum_per_snapshot = []

        # Load the simulation data
        data = sw.load(snapshot)

        boxsize = data.metadata.boxsize

        # Compute gas momentum
        gas_momentum = np.vstack((data.gas.masses, data.gas.masses, data.gas.masses)).T * data.gas.velocities

        pos = data.gas.coordinates - boxsize/2
        vel = data.gas.velocities
        r = np.linalg.norm(pos, axis=1)

        # Compute radial positions centered on explosion
        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

        # Compute radial velocity component
        p_r = np.abs((x * gas_momentum[:, 0] + y * gas_momentum[:, 1] + z * gas_momentum[:, 2]) / r)

        #%% Now plot quantities
        fig, ax = plt.subplots()
        ax.scatter(r, p_r.to(u.Msun * u.km/u.s), s=1, alpha=0.5)
        ax.set_xlabel(r"$r$ [kpc]")
        ax.set_ylabel(r"p$_r$ [M$_\odot$ km/s]")
        ax.grid(True, linestyle="--", alpha=0.6)

        # Save plot
        fig.tight_layout()
        plt.savefig(output_filename, format="png", bbox_inches='tight', dpi=300)
