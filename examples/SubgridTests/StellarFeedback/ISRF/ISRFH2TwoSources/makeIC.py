"""Initial conditions of two stars in a box with a dense molecular slab.

A slab of density contrast `--contrast` at a uniform H2 ratio shadows star A. Star B sits in
clear gas on the other side of the slab. `--stars` selects the stars written to the ICs: both
(the default), or one of them. The geometry file keeps both positions and the H2 densities,
which the check uses for the straight-ray column of each star.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import h2_fixture_common as common  # noqa: E402


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
        "--stars",
        choices=["ab", "a", "b"],
        default="ab",
        help="Stars written to the ICs",
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
    star_pos = {"ab": [star_a, star_b], "a": [star_a], "b": [star_b]}[o.stars]
    star_ids = {"ab": [n, n + 1], "a": [n], "b": [n + 1]}[o.stars]
    common.write_ic(
        o.out,
        box,
        pos,
        m0 * scale,
        rho0 * scale,
        np.zeros([n, 3]),
        np.array(star_pos),
        o.star_mass,
        star_ids=star_ids,
    )

    n_h = 0.76 * o.rho
    n_clear = o.h2_ratio * n_h
    common.write_fixture(
        o.fixture,
        kind="two_sources",
        box_kpc=box,
        star_a_kpc=star_a.tolist(),
        star_b_kpc=star_b.tolist(),
        slab_x_kpc=[x0, x1],
        n_h2_clear_cm3=n_clear,
        n_h2_slab_cm3=n_clear * o.contrast,
        n_gas=n,
    )
    print(f"Box {box:.4f} kpc, {n} gas, stars {o.stars}")


if __name__ == "__main__":
    main()
