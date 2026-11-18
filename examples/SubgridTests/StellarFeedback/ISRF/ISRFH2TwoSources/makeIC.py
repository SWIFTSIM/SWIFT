"""Initial conditions of two stars that light the same gas through different columns.

A dense molecular slab (density contrast at a uniform H2 ratio) shadows star A. Star B
sits in clear gas downstream of the slab. Gas close to star B receives both lights.
The gas is isothermal at 163.4 K, where the shielding fit does not depend on the density,
so the merge bounds apply. Set SPH:initial_temperature to the same value.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import h2_fixture_common as common  # noqa: E402


def min_image(d, box):
    return d - box * np.round(d / box)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--level", type=int, default=5)
    p.add_argument(
        "--rho", type=float, default=1e3, help="Clear-gas density in atom/cm3"
    )
    p.add_argument(
        "--mass", type=float, default=0.1, help="Clear-gas particle mass in Msun"
    )
    p.add_argument(
        "--contrast", type=float, default=300.0, help="Slab over clear density"
    )
    p.add_argument("--star_mass", type=float, default=29.7)
    p.add_argument("--h2_ratio", type=float, default=3e-4, help="nH2 / nH of the gas")
    p.add_argument(
        "--probe_radius_box",
        type=float,
        default=0.40,
        help="Probes lie within this distance of star B, in box sizes",
    )
    p.add_argument(
        "--merge_radius_box",
        type=float,
        default=0.12,
        help="The two lights merge at this distance from star B, in box sizes",
    )
    p.add_argument(
        "--temperature", type=float, default=163.4, help="Gas temperature in K"
    )
    p.add_argument("-o", dest="out", default="ICs_isrf_h2_two_sources.hdf5")
    p.add_argument("--fixture", default="fixture.json")
    o = p.parse_args()

    box, pos, n = common.box_positions(o.level, o.mass, o.rho)
    x0, x1 = 0.30 * box, 0.45 * box
    in_slab = (pos[:, 0] >= x0) & (pos[:, 0] < x1)
    m0 = o.mass * common.MSUN_CGS / common.UNIT_MASS_CGS
    rho0 = o.rho * common.M_PROTON_CGS * common.KPC_CGS**3 / common.UNIT_MASS_CGS
    scale = np.where(in_slab, o.contrast, 1.0)
    star_a = np.array([0.15 * box, 0.5 * box, 0.5 * box])
    star_b = np.array([0.60 * box, 0.5 * box, 0.5 * box])
    common.write_ic(
        o.out,
        box,
        pos,
        m0 * scale,
        rho0 * scale,
        np.zeros([n, 3]),
        np.array([star_a, star_b]),
        o.star_mass,
    )

    n_h = 0.76 * o.rho
    n_clear = o.h2_ratio * n_h
    n_slab = n_clear * o.contrast
    cm = common.KPC_CGS
    probes = []
    for i in range(n):
        db = min_image(pos[i] - star_b, box)
        rb = np.linalg.norm(db)
        if rb > o.probe_radius_box * box or rb == 0.0:
            continue
        da = pos[i] - star_a
        ra = np.linalg.norm(da)
        # Column of A's light: slab crossing along the ray plus clear gas for the rest.
        if da[0] > 0:
            cross = (
                max(0.0, min(da[0], x1 - star_a[0]) - max(0.0, x0 - star_a[0])) / da[0]
            )
        else:
            cross = 0.0
        path_slab = ra * cross
        n_a = (n_slab * path_slab + n_clear * (ra - path_slab)) * cm
        r_m = o.merge_radius_box * box
        down = n_clear * max(0.0, rb - r_m) * cm
        n_b = n_clear * min(rb, r_m) * cm
        probes.append(
            {
                "id": i,
                "r_a_kpc": ra,
                "r_b_kpc": rb,
                "downstream_N_cm2": down,
                "n_merges": 1,
                "shares": [
                    {"w": 1.0 / ra**2, "N_cm2": max(n_a - down, 0.0)},
                    {"w": 1.0 / min(rb, r_m) ** 2, "N_cm2": n_b},
                ],
            }
        )
    common.write_fixture(
        o.fixture,
        kind="two_sources",
        box_kpc=box,
        star_a_kpc=star_a.tolist(),
        star_b_kpc=star_b.tolist(),
        slab_x_kpc=[x0, x1],
        n_h2_clear_cm3=n_clear,
        n_h2_slab_cm3=n_slab,
        share_model="1/r^2 dilution, equal luminosities",
        temperature_K=o.temperature,
        merge_radius_kpc=o.merge_radius_box * box,
        probes=probes,
        n_gas=n,
    )
    print(f"Box {box:.4f} kpc, {n} gas, {len(probes)} probes near star B")


if __name__ == "__main__":
    main()
