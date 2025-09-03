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
# Parent directory containing multiple runs
base_folders = ["./"]

all_simulations_momentum = []

for base_folder in base_folders:
    # run_folders = sorted([os.path.join(base_folder, d) for d in os.listdir(
        # base_folder) if os.path.isdir(os.path.join(base_folder, d))])
    run_folders = ["./"]
    momentum_per_snapshot = []

# E_ej = 5.029144e-05, E_tot = 5.029144e-05, U_tot = 4.812343e-05, E_kin_tot = 1.408160e-05, p_ej = 2.350068e-07, p_terminal = 2.350068e-07, dU = 9.118633e-10, f_therm = 9.568910e-01
    p_ej = 2.3e-07*u.Msun * u.km / u.s
    
    # .to(u.Msun * u.km/u.s)

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
            
            momentum_per_snapshot = np.zeros((len(snapshots), 3))
            time_per_snapshot = np.zeros((len(snapshots), 1))
            E_kin_per_snapshot = np.zeros((len(snapshots), 1))
            E_int_per_snapshot = np.zeros((len(snapshots), 1))

            for index, snapshot in enumerate(snapshots):
                data = sw.load(snapshot)
                
                gas_momentum = np.vstack((data.gas.masses,data.gas.masses,data.gas.masses)).T * data.gas.velocities
                star_momentum = np.vstack((data.stars.masses,data.stars.masses,data.stars.masses)).T * data.stars.velocities
                momentum_tot = np.sum(gas_momentum + star_momentum, axis=0)

                E_kin = np.linalg.norm(gas_momentum)**2 / (2*data.gas.masses) + np.linalg.norm(star_momentum)**2 / (2*data.stars.masses)
                E_int = data.gas.internal_energies
                
                # Note: This example is without gravity => no E_pot_grav
                
                E_kin_per_snapshot[index] = np.sum(E_kin)
                E_int_per_snapshot[index] = np.sum(E_int)
                momentum_per_snapshot[index] = momentum_tot
                time_per_snapshot[index] = data.metadata.time

#%% Now plot quantities
            momentum_norm = np.linalg.norm(momentum_per_snapshot, axis=1)
            delta_p = np.diff(momentum_per_snapshot, axis=0)
            delta_p_norm = np.linalg.norm(delta_p, axis=1)
            
            fig, axes = plt.subplots(ncols=4, nrows=1, figsize=(16, 4))
   
            ax = axes[0]
            ax.plot(time_per_snapshot[1:], delta_p_norm)
            ax.set_yscale("log")
            ax.set_xlabel(r"$t$ [Gyr]")
            ax.set_ylabel("||p$_\mathrm{tot}$||")
            ax.grid(True, linestyle="--", alpha=0.6)
            
            ax = axes[1]
            ax.plot(time_per_snapshot, E_kin_per_snapshot)
            ax.set_yscale("log")
            ax.set_xlabel(r"$t$ [Gyr]")
            ax.set_ylabel("E$_\mathrm{kin}$")
            ax.grid(True, linestyle="--", alpha=0.6)
            
            ax = axes[2]
            ax.plot(time_per_snapshot, E_int_per_snapshot)
            ax.set_yscale("log")
            ax.set_xlabel(r"$t$ [Gyr]")
            ax.set_ylabel("E$_\mathrm{int}$")
            ax.grid(True, linestyle="--", alpha=0.6)

            ax = axes[3]
            ax.plot(time_per_snapshot, E_kin_per_snapshot+E_int_per_snapshot)
            ax.set_yscale("log")
            ax.set_xlabel(r"$t$ [Gyr]")
            ax.set_ylabel("E$_\mathrm{tot}$")
            ax.grid(True, linestyle="--", alpha=0.6)
            
            fig.subplots_adjust(left=0.06, right=0.985, top=0.97,
                                bottom=0.12, hspace=0.25, wspace=0)
            fig.tight_layout()
            plt.savefig("momentum_check.png",
                        format="png", bbox_inches='tight', dpi=300)
                