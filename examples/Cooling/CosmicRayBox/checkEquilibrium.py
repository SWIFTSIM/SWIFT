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
"""Compare the end state of the box with a one-zone run of pygeajay.

The bar is the 10% change per sub-step that geajay allows in the energy and
in the electrons: two runs that sub-cycle differently can end that far apart.
"""

import sys
from glob import glob

import h5py
import numpy as np
import pygeajay

rate = 2.0e-16  # s^-1, GeajayCooling:cosmic_ray_constant_rate
bar = 0.1
species = ["HI", "HII", "HeI", "HeII", "HeIII", "e"]


def load(name):
    with h5py.File(name, "r") as f:
        g = f["PartType0"]
        state = {k: g[k][:].astype(np.float64) for k in species}
        state["density"] = g["Densities"][:].astype(np.float64)
        state["energy"] = g["InternalEnergies"][:].astype(np.float64)
        u = f["Units"].attrs
        units = [float(u[k][0]) for k in ("Unit mass in cgs (U_M)", "Unit length in cgs (U_L)", "Unit time in cgs (U_t)")]
        time = float(f["Header"].attrs["Time"][0])
    return time, state, units


def fluid(chemistry, state, index):
    fc = pygeajay.FluidContainer(chemistry, 1)
    rho = state["density"][index]
    fc["density"][:] = rho
    for k in species:
        fc[k + "_density"][:] = state[k][index] * rho
    fc["internal_energy"][:] = state["energy"][index]
    return fc


def temperature_and_xe(fc):
    return float(fc.calculate_temperature()[0]), float(fc["e_density"][0] / (fc["HI_density"][0] + fc["HII_density"][0]))


snapshots = sorted(glob("cosmicRayBox_*.hdf5"))
t0, first, (um, ul, ut) = load(snapshots[0])
t1, last, _ = load(snapshots[-1])

units = pygeajay.CodeUnits(density_units=um / ul**3, length_units=ul, time_units=ut)
options = pygeajay.Options(chemistry_enabled=1, radiative_cooling=1, species_network=1, metal_line_cooling=0,
                           uv_background=0, case_b_recombination=1, hydrogen_mass_fraction=0.76,
                           adiabatic_index=5.0 / 3.0, max_substeps=10000, fail_on_max_substeps=1)
chemistry = options.to_chemistry_data()
chemistry.initialize(units)
cosmic_rays = pygeajay.CosmicRays(cosmic_ray_input="constant", cosmic_ray_constant_rate=float(np.float32(rate)))

one_zone = fluid(chemistry, first, 0)
for step in range(int(round((t1 - t0) / 1e-3))):
    one_zone.solve_chemistry_cosmic_ray_reduced(1e-3, cosmic_rays)
T_ref, xe_ref = temperature_and_xe(one_zone)

values = np.array([temperature_and_xe(fluid(chemistry, last, i)) for i in range(len(last["density"]))])
T_box, xe_box = np.median(values, axis=0)

dT, dxe = abs(T_box / T_ref - 1.0), abs(xe_box / xe_ref - 1.0)
print("SWIFT     : T = %.6g K, x_e = %.6g" % (T_box, xe_box))
print("pygeajay  : T = %.6g K, x_e = %.6g" % (T_ref, xe_ref))
ok = np.isfinite([dT, dxe]).all() and dT <= bar and dxe <= bar
print("relative differences %.2e and %.2e, bar %.2f: %s" % (dT, dxe, bar, "passed" if ok else "failed"))
sys.exit(0 if ok else 1)
