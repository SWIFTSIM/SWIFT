#!/usr/bin/env python3

"""
Exact solution for the Brio-Wu test.

Taken from plotSolution.jl.
"""

import numpy as np


def exact_Brio_Wu(time):

    tfac = 10.0 * time
    xmin = -1.0
    xmax = 1.0
    N = 14
    mu0_1 = 1.0 / np.sqrt(4.0 * np.pi)

    x = np.zeros(N)
    solution = np.zeros(
        N,
        dtype=[
            ("rho", np.float32),
            ("p", np.float32),
            ("u", np.float32),
            ("v", np.float32),
            ("w", np.float32),
            ("Bx", np.float32),
            ("By", np.float32),
            ("Bz", np.float32),
        ],
    )

    x[0] = xmin
    x[1] = -0.18 * tfac
    x[2] = -0.08 * tfac
    x[3:6] = -0.03 * tfac
    x[6] = -0.005 * tfac
    x[7:9] = 0.06 * tfac
    x[9:11] = 0.147 * tfac
    x[11] = 0.33 * tfac
    x[12] = 0.36 * tfac
    x[13] = xmax

    solution["rho"][0:2] = 1.0
    solution["rho"][2:4] = 0.67623
    solution["rho"][4] = 0.827
    solution["rho"][5] = 0.775
    solution["rho"][6:8] = 0.6962
    solution["rho"][8:10] = 0.2352
    solution["rho"][10:12] = 0.117
    solution["rho"][12:14] = 0.125

    solution["p"][0:2] = 1.0
    solution["p"][2:4] = 0.447
    solution["p"][4] = 0.727219
    solution["p"][5] = 0.6
    solution["p"][6:10] = 0.5160
    solution["p"][10:12] = 0.0876
    solution["p"][12:14] = 0.1

    solution["u"][0:2] = 0.0
    solution["u"][2:4] = 0.63721
    solution["u"][4] = 0.48
    solution["u"][5] = 0.52
    solution["u"][6:10] = 0.600
    solution["u"][10:12] = -0.24
    solution["u"][12:14] = 0.0

    solution["v"][0:2] = 0.0
    solution["v"][2:4] = -0.23345
    solution["v"][4] = -1.3
    solution["v"][5] = -1.4
    solution["v"][6:10] = -1.584
    solution["v"][10:12] = -0.166
    solution["v"][12:14] = 0.0

    solution["Bx"] = 0.75

    solution["By"][0:2] = 1.0
    solution["By"][2:4] = 2.1 * mu0_1
    solution["By"][4] = -1.2 * mu0_1
    solution["By"][5] = -1.3 * mu0_1
    solution["By"][6:10] = -1.9 * mu0_1
    solution["By"][10:12] = -3.25 * mu0_1
    solution["By"][12:14] = -1.0

    return x, solution
