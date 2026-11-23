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

This tests how SWIFT passes the particle state, the units and the cosmic-ray
rate to geajay. It does not test the cosmic-ray physics: SWIFT and pygeajay
call the same solver.

Usage: checkEquilibrium.py [snapshot_prefix] [--per-particle]

With --per-particle the rate is read from the snapshot field CosmicRayValues
and compared with the rate of the parameter file.

The end time is about 75 relaxation times after the gas reaches its
equilibrium, so the end state is the fixed point of the rate equations and does
not depend on the sub-stepping. SWIFT and pygeajay then differ only by the
float32 storage of the species and of the energy at every step (eps = 1.2e-7).
The bar is 1e-5, about 100 eps: it allows for the amplification of that
rounding by the temperature and the electron fraction, which is below 10. The
bar is applied to the largest deviation over the particles.

The stored rate is not amplified: CosmicRayValues is float32 and holds the
constant once rounded, so its bar is half an ulp, 2^-24 = 6e-8.

A wrong ionisation rate shifts the equilibrium: x_e changes by 0.5 d(rate) /
rate and T by 0.04 d(rate) / rate (pygeajay, rate +-1%, +-10%). The bar
detects an error of the rate above 2e-5 in x_e. It cannot detect an error that
SWIFT and pygeajay share, such as the heat per ionisation.
"""

import sys
from glob import glob

import h5py
import numpy as np
import pygeajay

rate = 2.0e-16  # s^-1, the rate of the parameter file
bar = 1.0e-5
bar_values = 2.0**-24  # half an ulp of float32: the stored rate is one rounding of its value
species = ["HI", "HII", "HeI", "HeII", "HeIII", "e"]


def load(name):
    with h5py.File(name, "r") as f:
        g = f["PartType0"]
        state = {k: g[k][:].astype(np.float64) for k in species}
        state["density"] = g["Densities"][:].astype(np.float64)
        state["energy"] = g["InternalEnergies"][:].astype(np.float64)
        if "CosmicRayValues" in g:
            state["cosmic_ray_value"] = g["CosmicRayValues"][:].astype(np.float64) * g["CosmicRayValues"].attrs["Conversion factor to CGS (not including cosmological corrections)"][0]
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


per_particle = "--per-particle" in sys.argv
arguments = [a for a in sys.argv[1:] if not a.startswith("--")]
prefix = arguments[0] if arguments else "cosmicRayBox"

snapshots = sorted(glob(prefix + "_*.hdf5"))
t0, first, (um, ul, ut) = load(snapshots[0])
t1, last, _ = load(snapshots[-1])

ok = True
if per_particle:
    values = last["cosmic_ray_value"]
    deviation = np.max(np.abs(values / rate - 1.0))
    print("CosmicRayValues: largest deviation from %g s^-1 is %.2e, bar %.1e" % (rate, deviation, bar_values))
    ok = ok and np.isfinite(deviation) and deviation <= bar_values

units = pygeajay.CodeUnits(density_units=um / ul**3, length_units=ul, time_units=ut)
options = pygeajay.Options(chemistry_enabled=1, radiative_cooling=1, species_network=1, metal_line_cooling=0,
                           uv_background=0, case_b_recombination=1, hydrogen_mass_fraction=0.76,
                           adiabatic_index=5.0 / 3.0, max_substeps=10000, fail_on_max_substeps=1)
chemistry = options.to_chemistry_data()
chemistry.initialize(units)
cosmic_rays = pygeajay.CosmicRays(cosmic_ray_input="constant", cosmic_ray_constant_rate=float(np.float32(rate)))


def one_zone_state(density_factor):
    """Run the one-zone reference to the end time at a scaled density."""
    fc = fluid(chemistry, first, 0)
    fc["density"][:] *= density_factor
    for k in species:
        fc[k + "_density"][:] *= density_factor
    for step in range(int(round((t1 - t0) / 1e-3))):
        fc.solve_chemistry_cosmic_ray_reduced(1e-3, cosmic_rays)
    return np.array(temperature_and_xe(fc))


# The SPH density of the lattice differs between particles by the float32
# rounding of its neighbour sum. The equilibrium depends on the density, so
# the first-order dependence is taken out with the slope of the reference.
step_ln_rho = 1.0e-3
reference = one_zone_state(1.0)
slope = (one_zone_state(1.0 + step_ln_rho) / reference - 1.0) / step_ln_rho
ln_rho = last["density"] / first["density"][0] - 1.0
print("largest density difference to the first particle %.2e, d ln(T, x_e) / d ln rho = %.3f, %.3f" % (np.max(np.abs(ln_rho)), slope[0], slope[1]))

values = np.array([temperature_and_xe(fluid(chemistry, last, i)) for i in range(len(last["density"]))])
deviation = np.abs(values / reference - 1.0 - np.outer(ln_rho, slope))
dT, dxe = np.max(deviation, axis=0)
print("SWIFT     : T = %.6g K, x_e = %.6g (median)" % tuple(np.median(values, axis=0)))
print("pygeajay  : T = %.6g K, x_e = %.6g" % tuple(reference))
ok = ok and np.isfinite([dT, dxe]).all() and dT <= bar and dxe <= bar
print("largest relative differences over the particles %.2e and %.2e, bar %.0e: %s" % (dT, dxe, bar, "passed" if ok else "failed"))
sys.exit(0 if ok else 1)
