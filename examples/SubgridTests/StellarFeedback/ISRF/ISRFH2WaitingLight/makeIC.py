"""Initial conditions of one star in a periodic box of molecular gas with dense clumps.

Each clump has a Gaussian density profile of peak contrast `--contrast` over the clear
gas and a 1/e radius of `--radius` box sizes. The clumps sit between `--r_min` and
`--r_max` box sizes from the star, which is at the box centre. The H2 ratio is uniform,
so the H2 density follows the gas density. The geometry file keeps the star and the clumps.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import h2_fixture_common as common  # noqa: E402


def place_clumps(n_clumps, radius, r_min, r_max, seed):
    """Return clump centres (box units) at least two radii apart from each other."""
    rng = np.random.default_rng(seed)
    centres = []
    while len(centres) < n_clumps:
        c = 0.5 + rng.uniform(-r_max, r_max, 3)
        r = np.linalg.norm(c - 0.5)
        if not r_min <= r <= r_max:
            continue
        if any(np.linalg.norm(c - o) < 2.0 * radius for o in centres):
            continue
        centres.append(c)
    return np.array(centres)


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
        "--contrast", type=float, default=100.0, help="Clump peak over clear density"
    )
    p.add_argument("--n_clumps", type=int, default=8)
    p.add_argument(
        "--radius", type=float, default=0.09, help="Clump 1/e radius in box sizes"
    )
    p.add_argument(
        "--r_min",
        type=float,
        default=0.12,
        help="Nearest clump centre to the star, box sizes",
    )
    p.add_argument(
        "--r_max",
        type=float,
        default=0.32,
        help="Farthest clump centre from the star, box sizes",
    )
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--star_mass", type=float, default=29.7)
    p.add_argument("--h2_ratio", type=float, default=3e-5, help="nH2 / nH of the gas")
    p.add_argument("-o", dest="out", default="ICs_isrf_h2_waiting_light.hdf5")
    p.add_argument("--fixture", default="fixture.json")
    o = p.parse_args()

    box, pos, n = common.box_positions(o.level, o.mass, o.rho)
    centres = place_clumps(o.n_clumps, o.radius, o.r_min, o.r_max, o.seed)
    scale = np.ones(n)
    for c in centres:
        d = pos / box - c[None, :]
        r2 = np.sum(d**2, axis=1)
        scale += (o.contrast - 1.0) * np.exp(-r2 / o.radius**2)
    m0 = o.mass * common.MSUN_CGS / common.UNIT_MASS_CGS
    rho0 = o.rho * common.M_PROTON_CGS * common.KPC_CGS**3 / common.UNIT_MASS_CGS
    star = np.full([1, 3], box / 2)
    common.write_ic(
        o.out, box, pos, m0 * scale, rho0 * scale, np.zeros([n, 3]), star, o.star_mass
    )
    n_h = 0.76 * o.rho
    common.write_fixture(
        o.fixture,
        kind="clumpy_cloud",
        box_kpc=box,
        star_kpc=star[0].tolist(),
        clump_centres_box=centres.tolist(),
        clump_radius_box=o.radius,
        contrast=o.contrast,
        n_h2_clear_cm3=o.h2_ratio * n_h,
        n_gas=n,
    )
    print(
        f"Box {box:.4f} kpc, {n} gas, {o.n_clumps} clumps, peak mass ratio {scale.max():.1f}"
    )


if __name__ == "__main__":
    main()
