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

# Checks that the Zeldovich pancake evolves adiabatically before shell
# crossing (caustic formation at z_c = 1).
#
# Before the caustic, the flow is a smooth compression in an expanding
# universe: there are no shocks and the entropy A = (gamma - 1) u / rho^(gamma-1)
# of every particle must be conserved. Any change is spurious heating (or
# cooling) by the artificial dissipation terms, e.g. viscosity acting on pairs
# that are separating in the physical (peculiar + Hubble) frame.
#
# Usage: python3 checkAdiabaticity.py [basename] [z_check] [tolerance]
# Returns a non-zero exit code if the entropy of any particle has changed by
# more than the tolerance by redshift z_check.

import glob
import sys
import h5py
import numpy as np

basename = sys.argv[1] if len(sys.argv) > 1 else "zeldovichPancake"
z_check = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
tolerance = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-2


def read(filename):
    with h5py.File(filename, "r") as f:
        gas = f["PartType0"]
        order = np.argsort(gas["ParticleIDs"][:])
        u = gas["InternalEnergies"][:][order]
        rho = gas["Densities"][:][order]
        gamma = float(np.atleast_1d(f["HydroScheme"].attrs["Adiabatic index"])[0])
        z = float(np.atleast_1d(f["Header"].attrs["Redshift"])[0])
    # Co-moving entropy (conserved along adiabatic trajectories)
    return z, (gamma - 1.0) * u / rho ** (gamma - 1.0)


files = sorted(glob.glob("%s_[0-9][0-9][0-9][0-9].hdf5" % basename))
if len(files) < 2:
    print("Need at least two snapshots named %s_XXXX.hdf5" % basename)
    sys.exit(1)

z_0, A_0 = read(files[0])
print("Entropy conservation before the caustic (reference: z = %.2f)" % z_0)
print("%10s %20s %20s" % ("z", "max |dA / A|", "median |dA / A|"))

worst = 0.0
checked = 0
for filename in files[1:]:
    z, A = read(filename)
    if z < z_check:
        break
    change = np.abs(A / A_0 - 1.0)
    worst = max(worst, change.max())
    checked += 1
    print("%10.3f %20.3e %20.3e" % (z, change.max(), np.median(change)))

if checked == 0:
    print("No snapshot with z >= %.2f to check" % z_check)
    sys.exit(1)

if worst > tolerance:
    print(
        "FAILED: the entropy changed by up to %.3e (tolerance %.1e) before z = %.2f"
        % (worst, tolerance, z_check)
    )
    sys.exit(1)

print(
    "PASSED: max entropy change %.3e < %.1e up to z = %.2f"
    % (worst, tolerance, z_check)
)
