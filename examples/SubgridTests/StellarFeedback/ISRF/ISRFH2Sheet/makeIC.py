"""Initial conditions of the H2 sheet: a square lattice of stars in the mid-plane of a uniform periodic box.

The box is two glass cubes stacked along z, so it is L x L x 2L. The stars sit in the plane
z = L at the centres of an n x n grid. The geometry file fixture.json holds the box, the
star positions and the gas density.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import h2_fixture_common as common  # noqa: E402

PC_KPC = 1e-3


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--level",
        type=int,
        default=5,
        help="(2**level)**3 gas particles per glass cube",
    )
    p.add_argument(
        "--rho",
        type=float,
        default=1316.0,
        help="Gas density in atom/cm3 (n_H = 0.76 rho)",
    )
    p.add_argument(
        "--mass", type=float, default=0.03825, help="Gas particle mass in Msun"
    )
    p.add_argument("--star_mass", type=float, default=4.9, help="Star mass in Msun")
    p.add_argument(
        "--star_spacing",
        type=float,
        default=0.84444,
        help="Star spacing in pc, rounded so that a whole number of stars fits the box side",
    )
    p.add_argument("-o", dest="out", default="ICs_isrf_h2_sheet.hdf5")
    p.add_argument("--fixture", default="fixture.json")
    o = p.parse_args()

    box, pos, n = common.box_positions(o.level, o.mass, o.rho)
    pos = np.concatenate([pos, pos + np.array([0.0, 0.0, box])])
    n_side = max(1, int(round(box / (o.star_spacing * PC_KPC))))
    a = box / n_side
    g = (np.arange(n_side) + 0.5) * a
    star = np.array([[x, y, box] for x in g for y in g])

    m = o.mass * common.MSUN_CGS / common.UNIT_MASS_CGS
    rho = o.rho * common.M_PROTON_CGS * common.KPC_CGS**3 / common.UNIT_MASS_CGS
    common.write_ic(
        o.out,
        [box, box, 2 * box],
        pos,
        np.full(2 * n, m),
        np.full(2 * n, rho),
        np.zeros([2 * n, 3]),
        star,
        o.star_mass,
        star_h_factor=1.0,
    )
    common.write_fixture(
        o.fixture,
        kind="sheet",
        box_kpc=[box, box, 2 * box],
        sheet_z_kpc=box,
        star_spacing_kpc=a,
        stars_per_side=n_side,
        star_kpc=star.tolist(),
        star_mass_msun=o.star_mass,
        n_hydrogen_cm3=0.76 * o.rho,
        n_gas=2 * n,
    )
    print(
        f"Box {box / PC_KPC:.3f} x {box / PC_KPC:.3f} x {2 * box / PC_KPC:.3f} pc, "
        f"{2 * n} gas particles, {len(star)} stars {a / PC_KPC:.4f} pc apart"
    )


if __name__ == "__main__":
    main()
