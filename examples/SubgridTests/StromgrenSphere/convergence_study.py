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
import os
import swiftsimio as sw
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
import h5py

# %%


def list_snapshots(folder_path):
    """Lists all snapshot files in the folder."""
    snapshots = [
        os.path.join(folder_path, f)
        for f in os.listdir(folder_path)
        if f.startswith("snapshot_") and f.endswith(".hdf5")
    ]
    return sorted(snapshots)




# %%
if __name__ == "__main__":
    # Parent directory containing multiple runs
    base_folders = ["./GEAR_MECH2"]
    plot_labels = ["Mechanical 2", "GEAR", "Mechanical 2"]
    colors = ["red", "blue", "orange"]
    simulation_names = ["Mechanical 2", "GEAR", "Mechanical 2"]
    data_output_filename = "isotropy_check_data.hdf5"
    # base_folders = ["./GEAR_MECH1", "./GEAR", "./GEAR_MECH2"]
    # plot_labels = ["Mechanical 1", "GEAR", "Mechanical 2"]
    # colors = ["red", "blue", "orange"]
    # simulation_names = ["Mechanical 1", "GEAR", "Mechanical 2"]
    # data_output_filename = "isotropy_check_data.hdf5"

    all_simulations_metal_flux_hists = []

    for base_folder in base_folders:
        run_folders = sorted([os.path.join(base_folder, d) for d in os.listdir(
            base_folder) if os.path.isdir(os.path.join(base_folder, d))])
        all_metal_flux_hists = []

        if not run_folders:
            print(f"No simulation folders found in {base_folder}.")
        else:
            print(f"Found {len(run_folders)} simulation runs.")

            for run_folder in tqdm(run_folders, desc="Processing simulations"):
                folder_path = os.path.join(run_folder, "snap")
                snapshots = list_snapshots(folder_path)

                if not snapshots:
                    print(f"No snapshots found in {folder_path}. Skipping.")
                    continue

                # Find first snapshot where feedback occurs
                for index, snapshot in enumerate(snapshots):
                    data = sw.load(snapshot)
                    if np.sum(data.gas.metal_mass_fractions.metals) > 0:

                        # Find the particles ids that received feedback at injection time and get the metals
                        I = np.argwhere(
                            data.gas.metal_mass_fractions.metals > 0).flatten()
                        m_metals = data.gas.metal_mass_fractions.metals[I] * \
                            data.gas.masses[I]
                        ids = data.gas.particle_ids[I]

                        # Get the snapshots
                        snapshot_injection = snapshot
                        snapshot_injection_1 = snapshots[index-1]
                        break

                ###################################
                # Now do data analysis

                # Load the data
                data = sw.load(snapshot_injection)

                

                # Do an histogram
                cos_theta_midpoints, metal_flux_hist = compute_flux_by_cos_theta(
                    star_pos, pos[mask], m_metals_flux[mask], 20)
                all_metal_flux_hists.append(metal_flux_hist)

            all_simulations_metal_flux_hists.append(all_metal_flux_hists)

# %%
    # Plot
    fig, ax = plt.subplots()
    ax.set_xlabel(r"$|\cos(\theta)|$", fontsize=14)
    ax.set_ylabel("Metal Flux", fontsize=14)
    for metal_flux_hist, plot_label, color in zip(all_simulations_metal_flux_hists, plot_labels, colors):

        # Compute mean and standard deviation
        median_flux_hist = np.median(metal_flux_hist, axis=0)
        std_flux_hist = np.std(metal_flux_hist, axis=0)

        # Plot
        ax.plot(cos_theta_midpoints, median_flux_hist,
                label=plot_label, color=color, linewidth=2)

    ax.legend(fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()
    plt.show()
    plt.savefig("isotropy_check_comparison.png",
                format="png", bbox_inches='tight', dpi=300)

    # Write the data
    with h5py.File(data_output_filename, "w") as f:
        for i, simulation in enumerate(simulation_names):
            simulation_group = f.create_group(simulation)

            if "metal_flux_hist" in simulation_group:
                # Assign values from output
                f.create_dataset("metal_flux_hist",
                                 data=all_simulations_metal_flux_hists[i])
                f.create_dataset("cos_theta", data=cos_theta_midpoints)