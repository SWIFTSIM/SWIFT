"""Initial conditions of the star-lit molecular box: one star in a uniform periodic box."""

import argparse

import numpy as np

import h2_fixture_common as common


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--level", type=int, default=5, help="N = (2**level)**3 gas particles"
    )
    p.add_argument("--rho", type=float, default=1e3, help="Density in atom/cm3")
    p.add_argument("--mass", type=float, default=0.1, help="Gas particle mass in Msun")
    p.add_argument("--star_mass", type=float, default=29.7, help="Star mass in Msun")
    p.add_argument("--h2_ratio", type=float, default=3e-4, help="nH2 / nH of the gas")
    p.add_argument("-o", dest="out", default="ICs_isrf_h2_star_lit.hdf5")
    p.add_argument("--fixture", default="fixture.json")
    o = p.parse_args()

    box, pos, n = common.box_positions(o.level, o.mass, o.rho)
    m = o.mass * common.MSUN_CGS / common.UNIT_MASS_CGS
    rho = o.rho * common.M_PROTON_CGS * common.KPC_CGS**3 / common.UNIT_MASS_CGS
    star = np.full([1, 3], box / 2)
    common.write_ic(
        o.out,
        box,
        pos,
        np.full(n, m),
        np.full(n, rho),
        np.zeros([n, 3]),
        star,
        o.star_mass,
    )
    n_h = 0.76 * o.rho
    common.write_fixture(
        o.fixture,
        kind="star_lit",
        box_kpc=box,
        star_kpc=star[0].tolist(),
        n_hydrogen_cm3=n_h,
        n_h2_cm3=o.h2_ratio * n_h,
        h2_ratio=o.h2_ratio,
        n_gas=n,
    )
    print(f"Box {box:.4f} kpc, {n} gas particles, n_H2 = {o.h2_ratio * n_h:.3g} cm^-3")


if __name__ == "__main__":
    main()
