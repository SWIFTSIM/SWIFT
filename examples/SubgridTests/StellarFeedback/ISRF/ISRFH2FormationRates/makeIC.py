"""Initial conditions of the dark formation box: uniform, periodic gas and no star."""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import h2_fixture_common as common  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--level", type=int, default=4, help="N = (2**level)**3 gas particles"
    )
    p.add_argument(
        "--rho",
        type=float,
        default=1e3 / 0.76,
        help="Density in atom/cm3 (n_H = 0.76 times this)",
    )
    p.add_argument("--mass", type=float, default=0.1, help="Gas particle mass in Msun")
    p.add_argument("-o", dest="out", default="ICs_isrf_h2_formation.hdf5")
    o = p.parse_args()

    box, pos, n = common.box_positions(o.level, o.mass, o.rho)
    m = o.mass * common.MSUN_CGS / common.UNIT_MASS_CGS
    rho = o.rho * common.M_PROTON_CGS * common.KPC_CGS**3 / common.UNIT_MASS_CGS
    common.write_ic(
        o.out,
        box,
        pos,
        np.full(n, m),
        np.full(n, rho),
        np.zeros([n, 3]),
        np.zeros([0, 3]),
        0.0,
    )
    print(
        f"Box {box:.4f} kpc, {n} gas particles, n_H = {0.76 * o.rho:.4g} cm^-3, no star"
    )


if __name__ == "__main__":
    main()
