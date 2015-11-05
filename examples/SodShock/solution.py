###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

import random
from numpy import *

# Generates the analytical  solution for the Sod shock test case
# The script works for a given left (x<0) and right (x>0) state and computes the solution at a later time t.
# The code writes five files rho.dat, P.dat, v.dat, u.dat and s.dat with the density, pressure, internal energy and
# entropic function on N points between x_min and x_max.
# This follows the solution given in (Toro, 2009)


# Parameters
rho_L = 1
P_L = 1
v_L = 0.

rho_R = 0.25
P_R = 0.1795
v_R = 0.

gamma = 5./3. # Polytropic index

t = 0.12  # Time of the evolution

N = 1000  # Number of points
x_min = -0.25
x_max = 0.25


# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

c_L = sqrt(gamma * P_L / rho_L)   # Speed of the rarefaction wave
c_R = sqrt(gamma * P_R / rho_R)   # Speed of the shock front

# Helpful variable
Gama = (gamma - 1.) / (gamma + 1.)
beta = (gamma - 1.) / (2. * gamma)

# Characteristic function and its derivative, following Toro (2009)
def compute_f(P_3, P, c):
    u = P_3 / P
    if u > 1:
        term1 = gamma*((gamma+1.)*u + gamma-1.)
        term2 = sqrt(2./term1)
        fp = (u - 1.)*c*term2
        dfdp = c*term2/P + (u - 1.)*c/term2*(-1./term1**2)*gamma*(gamma+1.)/P
    else:
        fp = (u**beta - 1.)*(2.*c/(gamma-1.))
        dfdp = 2.*c/(gamma-1.)*beta*u**(beta-1.)/P
    return (fp, dfdp)

# Solution of the Riemann problem following Toro (2009) 
def RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R):
    P_new = ((c_L + c_R + (v_L - v_R)*0.5*(gamma-1.))/(c_L / P_L**beta + c_R / P_R**beta))**(1./beta)
    P_3 = 0.5*(P_R + P_L)
    f_L = 1.
    while fabs(P_3 - P_new) > 1e-6:
        P_3 = P_new
        (f_L, dfdp_L) = compute_f(P_3, P_L, c_L)
        (f_R, dfdp_R) = compute_f(P_3, P_R, c_R)
        f = f_L + f_R + (v_R - v_L)
        df = dfdp_L + dfdp_R
        dp =  -f/df
        prnew = P_3 + dp
    v_3 = v_L - f_L
    return (P_new, v_3)


# Solve Riemann problem for post-shock region
(P_3, v_3) = RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R)

# Check direction of shocks and wave
shock_R = (P_3 > P_R)
shock_L = (P_3 > P_L)

# Velocity of shock front and and rarefaction wave
if shock_R:
    v_right = v_R + c_R**2*(P_3/P_R - 1.)/(gamma*(v_3-v_R))
else:
    v_right = c_R + 0.5*(gamma+1.)*v_3 - 0.5*(gamma-1.)*v_R

if shock_L:
    v_left = v_L + c_L**2*(P_3/p_L - 1.)/(gamma*(v_3-v_L))
else:
    v_left = c_L - 0.5*(gamma+1.)*v_3 + 0.5*(gamma-1.)*v_L

# Compute position of the transitions
x_23 = -fabs(v_left) * t
if shock_L :
    x_12 = -fabs(v_left) * t
else:
    x_12 = -(c_L - v_L) * t

x_34 = v_3 * t

x_45 = fabs(v_right) * t
if shock_R:
    x_56 = fabs(v_right) * t
else:
    x_56 = (c_R + v_R) * t 


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
            rho_s[i] = rho_L*(Gama + P_3/P_L)/(1. + Gama * P_3/P_L)
            P_s[i] = P_3
            v_s[i] = v_3
        else:
            rho_s[i] = rho_L*(Gama * (0. - x_s[i])/(c_L * t) + Gama * v_L/c_L + (1.-Gama))**(2./(gamma-1.))
            P_s[i] = P_L*(rho_s[i] / rho_L)**gamma
            v_s[i] = (1.-Gama)*(c_L -(0. - x_s[i]) / t) + Gama*v_L
    if x_s[i] >= x_23 and x_s[i] < x_34:
        if shock_L:
            rho_s[i] = rho_L*(Gama + P_3/P_L)/(1+Gama * P_3/p_L)
        else:
            rho_s[i] = rho_L*(P_3 / P_L)**(1./gamma)
        P_s[i] = P_3
        v_s[i] = v_3
    if x_s[i] >= x_34 and x_s[i] < x_45:
        if shock_R:
            rho_s[i] = rho_R*(Gama + P_3/P_R)/(1. + Gama * P_3/P_R)
        else:
            rho_s[i] = rho_R*(P_3 / P_R)**(1./gamma)
        P_s[i] = P_3
        v_s[i] = v_3
    if x_s[i] >= x_45 and x_s[i] < x_56:
        if shock_R:
            rho_s[i] = rho_R
            P_s[i] = P_R
            v_s[i] = v_R
        else:
            rho_s[i] = rho_R*(Gama*(x_s[i])/(c_R*t) - Gama*v_R/c_R + (1.-Gama))**(2./(gamma-1.))
            P_s[i] = p_R*(rho_s[i]/rho_R)**gamma
            v_s[i] = (1.-Gama)*(-c_R - (-x_s[i])/t) + Gama*v_R
    if x_s[i] >= x_56:
        rho_s[i] = rho_R
        P_s[i] = P_R
        v_s[i] = v_R


# Additional arrays
u_s = P_s / (rho_s * (gamma - 1.))  #internal energy
s_s = P_s / rho_s**gamma # entropic function
        
#---------------------------------------------------------------
# Print arrays

savetxt("rho.dat", column_stack((x_s, rho_s)))
savetxt("P.dat", column_stack((x_s, P_s)))
savetxt("v.dat", column_stack((x_s, v_s)))
savetxt("u.dat", column_stack((x_s, u_s)))
savetxt("s.dat", column_stack((x_s, s_s)))
