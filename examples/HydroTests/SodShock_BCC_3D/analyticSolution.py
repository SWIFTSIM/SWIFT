###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2019 Josh Borrow (joshua.boorrow@durham.ac.uk)
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
##############################################################################

import numpy as np
import sys

sys.path.append("../")

from riemannSolver import RiemannSolver


def analytic(
    time,  # Simulation t
    gas_gamma=5.0 / 3.0,  # Polytropic index
    rho_L=1.0,  # Density left state
    rho_R=0.125,  # Density right state
    v_L=0.0,  # Velocity left state
    v_R=0.0,  # Velocity right state
    P_L=1.0,  # Pressure left state
    P_R=0.1,  # Pressure right state
    N=1000,
    x_min=-1.0,
    x_max=1.0,
):

    solver = RiemannSolver(gas_gamma)
    delta_x = (x_max - x_min) / N
    x_s = np.arange(0.5 * x_min, 0.5 * x_max, delta_x)
    rho_s, v_s, P_s, _ = solver.solve(rho_L, v_L, P_L, rho_R, v_R, P_R, x_s / time)

    # Additional arrays
    u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy
    s_s = P_s / rho_s ** gas_gamma  # entropic function

    return dict(x=x_s + 1.0, v=v_s, rho=rho_s, P=P_s, u=u_s, S=s_s)
