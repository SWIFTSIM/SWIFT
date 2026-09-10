"""
Plot the cumulative halo mass function from the z=0 HBT-HERONS catalogue.

For more example analyses see https://hbt-herons.strw.leidenuniv.nl/examples/overview/
"""

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

hbt_output_path  = "./outputs/HBT-HERONS"
swift_output_path = "./outputs/SWIFT"

# Find the last available HBT-HERONS output.
snap_dirs = [d for d in os.listdir(hbt_output_path) if d.isdigit()]
if not snap_dirs:
    raise RuntimeError(f"No HBT-HERONS catalogues found in {hbt_output_path}")
last_snap = max(int(d) for d in snap_dirs)
snap_str  = f"{last_snap:03d}"

# Read all per-rank files for this snapshot and concatenate.
snap0 = f"{hbt_output_path}/{snap_str}/SubSnap_{snap_str}.0.hdf5"
with h5py.File(snap0, "r") as f:
    n_files      = f["NumberOfFiles"][0]
    mass_unit    = f["Units/MassInMsunh"][0]   # internal mass unit in M_sun/h
    swift_snap   = f["SnapshotId"][0]

data = np.concatenate([
    h5py.File(f"{hbt_output_path}/{snap_str}/SubSnap_{snap_str}.{i}.hdf5", "r")["Subhalos"][()]
    for i in range(n_files)
])

# Read h from the corresponding SWIFT snapshot.
swift_snap_file = f"{swift_output_path}/snap_{swift_snap:04d}.hdf5"
with h5py.File(swift_snap_file, "r") as f:
    h = f["Cosmology"].attrs["h"]

# Convert mass from M_sun/h to M_sun.
Mbound = data["Mbound"] * mass_unit / h

# Select resolved central subhaloes (most massive subhalo per FoF group).
central_mask     = (data["Rank"] == 0) & (data["Nbound"] > 0)
Mbound_centrals  = np.sort(Mbound[central_mask])[::-1]

fig, ax = plt.subplots()
ax.plot(Mbound_centrals, np.arange(1, len(Mbound_centrals) + 1))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$M_{\rm bound} \; [\mathrm{M}_{\odot}]$")
ax.set_ylabel(r"$N(\geq M_{\rm bound})$")
ax.set_title(f"Halo mass function at snapshot {last_snap} (z $\\approx$ 0)")
fig.tight_layout()
fig.savefig("HMF.png", dpi=150)
print("Saved HMF.png")
