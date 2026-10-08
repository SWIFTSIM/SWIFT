################################################################################
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
################################################################################
"""Print the time-integration overrides for one physical duration.

Without cosmology (``--redshift 0``) the overrides are in internal time. With
cosmology, SWIFT integrates in ln a: ``a_end`` is the scale factor reached
after the requested proper time, ``dt_max`` is a step in ln a and the snapshot
``delta_time`` is a ratio of scale factors. The Friedmann equation uses the
``Cosmology`` block of ``params.yml`` (flat, Omega_r = 0), as SWIFT does.
The output is a list of shell assignments for ``run.sh``.
"""

import argparse

import numpy as np
import yaml
from scipy.integrate import quad
from scipy.optimize import brentq

KPC_PER_MPC = 1.0e3


def hubble_rate(a: float, cosmology: dict) -> float:
    """Return H(a) in internal units (km/s/kpc)."""
    omega_m = cosmology["Omega_cdm"] + cosmology["Omega_b"]
    omega_l = cosmology["Omega_lambda"]
    omega_r = cosmology.get("Omega_r", 0.0)
    omega_k = 1.0 - omega_m - omega_l - omega_r
    h0 = 100.0 * cosmology["h"] / KPC_PER_MPC
    return h0 * np.sqrt(omega_r * a**-4 + omega_m * a**-3 + omega_k * a**-2 + omega_l)


def proper_time(a_begin: float, a_end: float, cosmology: dict) -> float:
    """Return the proper time between two scale factors, internal units."""
    if a_end == a_begin:
        return 0.0
    # Integrate in ln a: the span can be as small as 1e-8.
    return quad(
        lambda x: 1.0 / hubble_rate(a_begin * np.exp(x), cosmology),
        0.0,
        np.log(a_end / a_begin),
        epsabs=0.0,
        epsrel=1e-13,
    )[0]


def main() -> None:
    """Print the overrides."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", default="params.yml")
    parser.add_argument("--redshift", type=float, required=True)
    parser.add_argument("--duration", type=float, required=True)
    parser.add_argument("--snapshots", type=int, required=True)
    parser.add_argument("--steps", type=int, required=True)
    opt = parser.parse_args()

    if opt.redshift <= 0.0:
        print(f"time_end={opt.duration:.17g}")
        print(f"dt_max={opt.duration / opt.steps:.17g}")
        print(f"delta_time={opt.duration / opt.snapshots:.17g}")
        return

    with open(opt.params) as handle:
        cosmology = yaml.safe_load(handle)["Cosmology"]
    a_begin = 1.0 / (1.0 + opt.redshift)
    h_begin = hubble_rate(a_begin, cosmology)
    upper = opt.duration * h_begin * 2.0 + 1e-12
    dlna = brentq(
        lambda x: proper_time(a_begin, a_begin * np.exp(x), cosmology) - opt.duration,
        0.0,
        upper,
        xtol=1e-16,
        rtol=1e-15,
    )
    a_end = a_begin * np.exp(dlna)
    print(f"a_begin={a_begin:.17g}")
    print(f"a_end={a_end:.17g}")
    print(f"dt_max={dlna / opt.steps:.17g}")
    print(f"delta_time={np.exp(dlna / opt.snapshots):.17g}")


if __name__ == "__main__":
    main()
