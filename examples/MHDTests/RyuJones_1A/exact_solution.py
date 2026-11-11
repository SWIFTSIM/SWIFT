#!/usr/bin/env python3

"""
Exact solution for the Ryu&Jones 1A shock tube test.
"""

import numpy as np


def exact_RyuJones_1A(time):

    tfac = 12.5 * time
    xmin = -1.0
    xmax = 1.0
    N = 12

    x = np.zeros(N)
    solution = np.zeros(
        N,
        dtype=[
            ("rho", np.float32),
            ("P", np.float32),
            # ("u_int", np.float32),
            ("u", np.float32),
            ("v", np.float32),
            ("w", np.float32),
            ("Bx", np.float32),
            ("By", np.float32),
        ],
    )

    x[0] = xmin
    x[1:3] = -0.386 * tfac
    x[3:5] = -0.01 * tfac
    x[5:7] = 0.0505 * tfac
    x[7:9] = 0.12 * tfac
    x[9:11] = 0.37 * tfac
    x[11] = xmax

    solution["rho"][0:2] = 1.0
    solution["rho"][2:4] = 2.6797
    solution["rho"][4:6] = 2.6713
    solution["rho"][6:8] = 3.8508
    solution["rho"][8:10] = 3.7481
    solution["rho"][10:12] = 1.0

    """
    solution["u_int"][0:2] = 30.0
    solution["u_int"][2:6] = 85.0
    solution["u_int"][6:8] = 84.0
    solution["u_int"][8:10] = 59.0
    solution["u_int"][10:12] = 58.0
    solution["u_int"][12:14] = 1.5
    """

    solution["P"][0:2] = 20.0
    solution["P"][2:4] = 150.98
    solution["P"][4:8] = 150.19
    solution["P"][8:10] = 143.57
    solution["P"][10:12] = 1.0

    solution["u"][0:2] = 10.0
    solution["u"][2:4] = 0.72113
    solution["u"][4:8] = 0.72376
    solution["u"][8:10] = 0.70505
    solution["u"][10:12] = -10.0

    solution["v"][0:2] = 0.0
    solution["v"][2:4] = 0.23139
    solution["v"][4:8] = 0.35684
    solution["v"][8:10] = -0.38804
    solution["v"][10:12] = 0.0

    B0 = 2.5 / np.sqrt(np.pi)

    solution["Bx"] = B0

    solution["By"][0:2] = B0
    solution["By"][2:4] = 3.8389
    solution["By"][4:8] = 4.0380
    solution["By"][8:10] = 5.4272
    solution["By"][10:12] = B0

    return x, solution
