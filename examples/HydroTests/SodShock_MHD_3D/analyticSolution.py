###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

from numpy import *


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
    # Analytic solution
    c_L = sqrt(gas_gamma * P_L / rho_L)  # Speed of the rarefaction wave
    c_R = sqrt(gas_gamma * P_R / rho_R)  # Speed of the shock front

    # Helpful variable
    Gama = (gas_gamma - 1.0) / (gas_gamma + 1.0)
    beta = (gas_gamma - 1.0) / (2.0 * gas_gamma)

    # Characteristic function and its derivative, following Toro (2009)
    def compute_f(P_3, P, c):
        u = P_3 / P
        if u > 1:
            term1 = gas_gamma * ((gas_gamma + 1.0) * u + gas_gamma - 1.0)
            term2 = sqrt(2.0 / term1)
            fp = (u - 1.0) * c * term2
            dfdp = (
                c * term2 / P
                + (u - 1.0)
                * c
                / term2
                * (-1.0 / term1 ** 2)
                * gas_gamma
                * (gas_gamma + 1.0)
                / P
            )
        else:
            fp = (u ** beta - 1.0) * (2.0 * c / (gas_gamma - 1.0))
            dfdp = 2.0 * c / (gas_gamma - 1.0) * beta * u ** (beta - 1.0) / P
        return (fp, dfdp)

    # Solution of the Riemann problem following Toro (2009)
    def RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R):
        P_new = (
            (c_L + c_R + (v_L - v_R) * 0.5 * (gas_gamma - 1.0))
            / (c_L / P_L ** beta + c_R / P_R ** beta)
        ) ** (1.0 / beta)
        P_3 = 0.5 * (P_R + P_L)
        f_L = 1.0
        while fabs(P_3 - P_new) > 1e-6:
            P_3 = P_new
            (f_L, dfdp_L) = compute_f(P_3, P_L, c_L)
            (f_R, dfdp_R) = compute_f(P_3, P_R, c_R)
            f = f_L + f_R + (v_R - v_L)
            df = dfdp_L + dfdp_R
            dp = -f / df
            prnew = P_3 + dp
        v_3 = v_L - f_L
        return (P_new, v_3)

    # Solve Riemann problem for post-shock region
    (P_3, v_3) = RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R)

    # Check direction of shocks and wave
    shock_R = P_3 > P_R
    shock_L = P_3 > P_L

    # Velocity of shock front and and rarefaction wave
    if shock_R:
        v_right = v_R + c_R ** 2 * (P_3 / P_R - 1.0) / (gas_gamma * (v_3 - v_R))
    else:
        v_right = c_R + 0.5 * (gas_gamma + 1.0) * v_3 - 0.5 * (gas_gamma - 1.0) * v_R

    if shock_L:
        v_left = v_L + c_L ** 2 * (P_3 / p_L - 1.0) / (gas_gamma * (v_3 - v_L))
    else:
        v_left = c_L - 0.5 * (gas_gamma + 1.0) * v_3 + 0.5 * (gas_gamma - 1.0) * v_L

    # Compute position of the transitions
    x_23 = -fabs(v_left) * time
    if shock_L:
        x_12 = -fabs(v_left) * time
    else:
        x_12 = -(c_L - v_L) * time

    x_34 = v_3 * time

    x_45 = fabs(v_right) * time
    if shock_R:
        x_56 = fabs(v_right) * time
    else:
        x_56 = (c_R + v_R) * time

    # Prepare arrays
    delta_x = (x_max - x_min) / N
    x_s = arange(x_min, x_max, delta_x)
    rho_s = zeros(N)
    P_s = zeros(N)
    v_s = zeros(N)

    # Compute solution in the different regions
    for i in range(N):
        if x_s[i] <= x_12:
            rho_s[i] = rho_L
            P_s[i] = P_L
            v_s[i] = v_L
        if x_s[i] >= x_12 and x_s[i] < x_23:
            if shock_L:
                rho_s[i] = rho_L * (Gama + P_3 / P_L) / (1.0 + Gama * P_3 / P_L)
                P_s[i] = P_3
                v_s[i] = v_3
            else:
                rho_s[i] = rho_L * (
                    Gama * (0.0 - x_s[i]) / (c_L * time)
                    + Gama * v_L / c_L
                    + (1.0 - Gama)
                ) ** (2.0 / (gas_gamma - 1.0))
                P_s[i] = P_L * (rho_s[i] / rho_L) ** gas_gamma
                v_s[i] = (1.0 - Gama) * (c_L - (0.0 - x_s[i]) / time) + Gama * v_L
        if x_s[i] >= x_23 and x_s[i] < x_34:
            if shock_L:
                rho_s[i] = rho_L * (Gama + P_3 / P_L) / (1 + Gama * P_3 / p_L)
            else:
                rho_s[i] = rho_L * (P_3 / P_L) ** (1.0 / gas_gamma)
            P_s[i] = P_3
            v_s[i] = v_3
        if x_s[i] >= x_34 and x_s[i] < x_45:
            if shock_R:
                rho_s[i] = rho_R * (Gama + P_3 / P_R) / (1.0 + Gama * P_3 / P_R)
            else:
                rho_s[i] = rho_R * (P_3 / P_R) ** (1.0 / gas_gamma)
            P_s[i] = P_3
            v_s[i] = v_3
        if x_s[i] >= x_45 and x_s[i] < x_56:
            if shock_R:
                rho_s[i] = rho_R
                P_s[i] = P_R
                v_s[i] = v_R
            else:
                rho_s[i] = rho_R * (
                    Gama * (x_s[i]) / (c_R * time) - Gama * v_R / c_R + (1.0 - Gama)
                ) ** (2.0 / (gas_gamma - 1.0))
                P_s[i] = p_R * (rho_s[i] / rho_R) ** gas_gamma
                v_s[i] = (1.0 - Gama) * (-c_R - (-x_s[i]) / time) + Gama * v_R
        if x_s[i] >= x_56:
            rho_s[i] = rho_R
            P_s[i] = P_R
            v_s[i] = v_R

    # Additional arrays
    u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy
    s_s = P_s / rho_s ** gas_gamma  # entropic function

    return dict(x=x_s + 1, v=v_s, rho=rho_s, P=P_s, u=u_s, S=s_s)
