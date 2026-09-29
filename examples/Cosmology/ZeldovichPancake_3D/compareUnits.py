################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
################################################################################

# Checks that the Zeldovich pancake does not depend on the internal unit
# system. It compares the run with zeldovichPancake.yml (Mpc) with the one with
# zeldovichPancake_kpc.yml (kpc), which describe the same physical problem.
#
# All quantities are converted to physical cgs units. Before the caustic,
# particles are compared one by one. After shell crossing, tiny differences are
# amplified by the collapse, so binned profiles along x are compared instead.
#
# Any dependence on the units reveals a hidden dimensional constant in the
# scheme (e.g. an absolute epsilon or a missing factor of h).
#
# Usage: python3 compareUnits.py [basename_1] [basename_2] [z_caustic_margin]

import glob
import sys
import h5py
import numpy as np

base_1 = sys.argv[1] if len(sys.argv) > 1 else "zeldovichPancake"
base_2 = sys.argv[2] if len(sys.argv) > 2 else "zeldovichPancake_kpc"
z_particle = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0

tol_particle = 1e-3  # per particle, before the caustic
tol_profile = 5e-2  # binned profiles (relative to the profile's range)
num_bins = 64

conversion = "Conversion factor to physical CGS (including cosmological corrections)"


def read(filename):
    with h5py.File(filename, "r") as f:
        gas = f["PartType0"]
        order = np.argsort(gas["ParticleIDs"][:])
        out = {"z": float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])}
        box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        for name in ("Coordinates", "Velocities", "Densities", "InternalEnergies"):
            d = gas[name]
            out[name] = d[:][order] * float(np.atleast_1d(d.attrs[conversion])[0])
        # Position along x as a fraction of the box
        x_conv = float(np.atleast_1d(gas["Coordinates"].attrs[conversion])[0])
        out["x"] = out["Coordinates"][:, 0] / (box * x_conv)
    return out


def profile(x, y):
    bins = np.linspace(0.0, 1.0, num_bins + 1)
    idx = np.clip(np.digitize(x, bins) - 1, 0, num_bins - 1)
    return np.array([np.mean(y[idx == b]) if np.any(idx == b) else np.nan for b in range(num_bins)])


files_1 = sorted(glob.glob("%s_[0-9][0-9][0-9][0-9].hdf5" % base_1))
files_2 = sorted(glob.glob("%s_[0-9][0-9][0-9][0-9].hdf5" % base_2))
n = min(len(files_1), len(files_2))
if n == 0:
    print("No snapshots to compare")
    sys.exit(1)

fields = [("Densities", "rho"), ("InternalEnergies", "u"), ("Velocities", "v_x")]
print("%8s %8s   %s" % ("snap", "z", "   ".join("%-24s" % f[1] for f in fields)))
failed = False
for i in range(n):
    s1, s2 = read(files_1[i]), read(files_2[i])
    per_particle = s1["z"] >= z_particle
    row = []
    for name, label in fields:
        y1 = s1[name][:, 0] if name == "Velocities" else s1[name]
        y2 = s2[name][:, 0] if name == "Velocities" else s2[name]
        if per_particle:
            scale = np.maximum(np.abs(y1), 1e-3 * np.abs(y1).max())
            err, tol, kind = np.max(np.abs(y2 - y1) / scale), tol_particle, "part."
        else:
            p1, p2 = profile(s1["x"], y1), profile(s2["x"], y2)
            err = np.nanmax(np.abs(p2 - p1)) / (np.nanmax(p1) - np.nanmin(p1))
            tol, kind = tol_profile, "prof."
        failed |= err > tol
        row.append("%s %.2e %-8s" % (kind, err, "" if err <= tol else "FAIL"))
    print("%8d %8.3f   %s" % (i, s1["z"], "   ".join(row)))

if failed:
    print("FAILED: the results depend on the internal unit system")
    sys.exit(1)
print("PASSED: the results do not depend on the internal unit system")
