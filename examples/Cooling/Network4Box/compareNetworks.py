###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
###############################################################################
"""Print the relative difference of network 4 against network 3 at the end of the box runs.

Usage: compareNetworks.py [prefix_of_network_3 [prefix_of_network_4]]

A missing prefix defaults to network3 or network4. The script takes the last
snapshot of each run and prints, for the median particle, the H2 fraction
(the share of the hydrogen mass in H2), the HD to H number ratio and the
temperature, with (network 4 / network 3 - 1). The temperature comes from
pygeajay, which needs GEAJAY_LIBRARY_PATH.
"""

import sys
from glob import glob

import h5py
import numpy as np
import pygeajay

x_h = 0.76
species = ["HI", "HII", "HeI", "HeII", "HeIII", "e", "HM", "H2I", "H2II", "DI", "DII", "HDI", "DM", "HDII", "HeHII"]


def load(prefix):
    """Return the median state of the last snapshot of a run, and its units."""
    name = sorted(glob(prefix + "_[0-9][0-9][0-9][0-9].hdf5"))[-1]
    with h5py.File(name, "r") as f:
        g = f["PartType0"]
        state = {k: float(np.median(g[k][:].astype(np.float64))) for k in species if k in g}
        state["density"] = float(np.median(g["Densities"][:].astype(np.float64)))
        state["energy"] = float(np.median(g["InternalEnergies"][:].astype(np.float64)))
        u = f["Units"].attrs
        units = [float(u[k][0]) for k in ("Unit mass in cgs (U_M)", "Unit length in cgs (U_L)", "Unit time in cgs (U_t)")]
    return state, units


def observables(state, units, network):
    """Return the H2 fraction, n(HD) / n_H and T of a median state."""
    um, ul, ut = units
    code_units = pygeajay.CodeUnits(density_units=um / ul**3, length_units=ul, time_units=ut)
    options = pygeajay.Options(chemistry_enabled=1, radiative_cooling=1, species_network=network, metal_line_cooling=0,
                               uv_background=0, case_b_recombination=1, hydrogen_mass_fraction=x_h,
                               adiabatic_index=5.0 / 3.0)
    chemistry = options.to_chemistry_data()
    chemistry.initialize(code_units)
    fc = pygeajay.FluidContainer(chemistry, 1)
    rho = state["density"]
    fc["density"][:] = rho
    for k in species:
        if k in state:
            fc[k + "_density"][:] = state[k] * rho
    fc["internal_energy"][:] = state["energy"]
    return state["H2I"] / x_h, state["HDI"] / 3.0 / x_h, float(fc.calculate_temperature()[0])


def compare(prefixes):
    """Return the observables of the networks 3 and 4 runs with the given prefixes."""
    rows = []
    for prefix, network in zip(prefixes, (3, 4)):
        state, units = load(prefix)
        rows.append(observables(state, units, network))
    return rows


if __name__ == "__main__":
    prefixes = sys.argv[1:3]
    rows = compare(prefixes + ["network3", "network4"][len(prefixes):])
    names = ["H2 fraction", "HD / H number ratio", "temperature [K]"]
    print("%-22s %14s %14s %12s" % ("", "network 3", "network 4", "difference"))
    for i, name in enumerate(names):
        print("%-22s %14.6g %14.6g %+11.3f%%" % (name, rows[0][i], rows[1][i], 100.0 * (rows[1][i] / rows[0][i] - 1.0)))
