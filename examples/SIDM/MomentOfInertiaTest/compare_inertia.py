###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
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
#################################################################################

# Compares the inertia tensors (and smoothing lengths) of SIDM and gas
# particles, matched on particle ID.

import argparse
import h5py
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--sidm", type=str, default="snap_sidm/snapshot_0000.hdf5")
parser.add_argument("--gas", type=str, default="snap_gas/snapshot_0000.hdf5")
args = parser.parse_args()


def read(fname, group):
    with h5py.File(fname, "r") as f:
        ids = f[group]["ParticleIDs"][:]
        order = np.argsort(ids)
        return (
            ids[order],
            f[group]["SmoothingLengths"][:][order],
            f[group]["InertiaTensors"][:][order],
        )

id_s, h_s, I_s = read(args.sidm, "SIDMParticles")
id_g, h_g, I_g = read(args.gas, "GasParticles")

print(f"Comparing {args.sidm} and {args.gas}")

dh = np.abs(h_s / h_g - 1.0)
print(f"Smoothing lengths: max |h_sidm/h_gas - 1| = {dh.max():.3e}")

# Normalise the differences by the trace of the tensor
trace = I_g[:, :3].sum(axis=1)
diff = np.abs(I_s - I_g) / trace[:, None]
labels = ["xx", "yy", "zz", "xy", "xz", "yz"]
for k, lab in enumerate(labels):
    print(f"I_{lab}: max |dI|/tr(I) = {diff[:, k].max():.3e}")

print(
    "Tensors approximately equal:",
    np.allclose(I_s, I_g, rtol=0.0, atol=1e-5 * trace[:, None]),
)
print("Number of particles with max |dI|/tr(I) > 1e-5:", (diff.max(axis=1) > 1e-5).sum())
