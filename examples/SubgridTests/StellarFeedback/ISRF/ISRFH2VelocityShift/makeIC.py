"""Initial conditions of the velocity-shift fixture.

A star lights the gas at rest inside a sphere around it, then a molecular cloud that starts
at a spherical face centred on the star and moves radially away from it. The cloud speed
relative to the inner gas is set in units of the thermal line width b of the H2 gas.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import h2_fixture_common as common  # noqa: E402

K_B = 1.380649e-16
M_H = 1.6735575e-24


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level", type=int, default=5)
    p.add_argument("--rho", type=float, default=1e3, help="Density in atom/cm3")
    p.add_argument("--mass", type=float, default=0.1)
    p.add_argument("--star_mass", type=float, default=29.7)
    p.add_argument("--h2_ratio", type=float, default=3e-4)
    p.add_argument(
        "--temperature", type=float, default=50.0, help="Gas temperature in K"
    )
    p.add_argument(
        "--widths",
        type=float,
        default=10.0,
        help="Cloud speed relative to the slab, in line widths b",
    )
    p.add_argument("-o", dest="out", default="ICs_isrf_h2_velocity_shift.hdf5")
    p.add_argument("--fixture", default="fixture.json")
    o = p.parse_args()

    box, pos, n = common.box_positions(o.level, o.mass, o.rho)
    star = np.array([0.15 * box, 0.5 * box, 0.5 * box])
    cloud = (0.25 * box, 0.90 * box)
    d = pos - star[None, :]
    d -= box * np.round(d / box)
    r = np.linalg.norm(d, axis=1)
    in_cloud = (r >= cloud[0]) & (r < cloud[1])
    b_kms = np.sqrt(2 * K_B * o.temperature / (2 * M_H)) / 1e5
    vel = np.zeros([n, 3])
    vel[in_cloud] = o.widths * b_kms * d[in_cloud] / r[in_cloud, None]
    m0 = o.mass * common.MSUN_CGS / common.UNIT_MASS_CGS
    rho0 = o.rho * common.M_PROTON_CGS * common.KPC_CGS**3 / common.UNIT_MASS_CGS
    common.write_ic(
        o.out,
        box,
        pos,
        np.full(n, m0),
        np.full(n, rho0),
        vel,
        star[None, :],
        o.star_mass,
    )
    n_h2 = o.h2_ratio * 0.76 * o.rho
    common.write_fixture(
        o.fixture,
        kind="velocity_shift",
        box_kpc=box,
        star_kpc=star.tolist(),
        cloud_radius_kpc=list(cloud),
        n_h2_cm3=n_h2,
        line_width_kms=b_kms,
        cloud_speed_widths=o.widths,
        temperature_K=o.temperature,
        n_gas=n,
    )
    print(f"Box {box:.4f} kpc, b = {b_kms:.3f} km/s, cloud speed {o.widths} b")


if __name__ == "__main__":
    main()
