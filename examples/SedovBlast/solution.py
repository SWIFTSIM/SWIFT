###############################################################################
 # This file is part of SWIFT.
 # Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

# Computes the analytical solution of the 3D Sedov blast wave.
# The script works for a given initial box and dumped energy and computes the solution at a later time t.
# The code writes five files rho.dat, P.dat, v.dat, u.dat and s.dat with the density, pressure, internal energy and
# entropic function on N radial points between r=0 and r=R_max.
# Follows the solution in Landau & Lifschitz

# Parameters
rho_0 = 1.          # Background Density
P_0 = 1.e-5         # Background Pressure
E_0 = 1.e2          # Energy of the explosion
gamma = 5./3.     # Gas polytropic index

t = 0.15           # Time of the solution

N = 1000          # Number of radial points
R_max = 3.        # Maximal radius

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

if gamma == 5./3.:
    alpha = 0.49 
else:
    print "Unknown value for alpha"
    exit

# Position and velocity of the shock
r_shock = (E_0  / (alpha * rho_0))**(1./5.) * t**(2./5.)
v_shock = (1./5.) * (1./alpha)**(1./5.) * ((E_0 * t**2. / rho_0)**(-4./5.)) * E_0 * (t / rho_0)

# Prepare arrays
delta_r = R_max / N
r_s = arange(0, R_max, delta_r) 
rho_s = ones(N) * rho_0
P_s = ones(N) * P_0
u_s = ones(N)
v_s = zeros(N)

# State on the shock
rho_shock = rho_0 * (gamma+1.)/(gamma-1.)
P_shock = 2./(gamma+1.) * rho_shock * v_shock**2

# Integer position of the shock 
i_shock = min(floor(r_shock /delta_r), N)

# Dimensionless velocity and its spatial derivative
v_bar0 = (gamma+1.)/(2.*gamma)
deltaV_bar = (1.0 - v_bar0) / (i_shock - 1)

def rho_dimensionless(v_bar):
    power1 = (2./(gamma-2.))
    power2 = -(12.-7.*gamma+13.*gamma**2)/(2.-3.*gamma-11.*gamma**2+6.*gamma**3)
    power3 = 3./(1.+2.*gamma)
    term1 = ((1.+ gamma - 2.*v_bar)/(gamma-1.))**power1
    term2 = ((5.+5.*gamma+2.*v_bar-6.*gamma*v_bar)/(7.-gamma))**power2
    term3 = ((2.*gamma*v_bar - gamma -1.)/(gamma-1.))**power3
    return term1 * term2 * term3

def P_dimensionless(v_bar):
    return (gamma+1. - 2.*v_bar)/(2.*gamma*v_bar - gamma - 1.)*v_bar**2 * rho_dimensionless(v_bar)

def r_dimensionless(v_bar):
    power1 = (-12.+7.*gamma-13.*gamma**2)/(5.*(-1.+gamma+6.*gamma**2))
    power2 = (gamma - 1.)/(1.+2.*gamma)
    term1 = ((5.+5.*gamma+2.*v_bar-6.*gamma*v_bar)/(7.-gamma))**power1
    term2 = ((2.*gamma*v_bar-gamma-1.)/(gamma-1.))**power2
    return v_bar**(-2./5.)*term1*term2

# Generate solution
for i in range(1,int(i_shock)):
    v_bar = v_bar0 + (i-1)*deltaV_bar    
    r_s[i] = r_dimensionless(v_bar) * r_shock
    rho_s[i] = rho_dimensionless(v_bar) * rho_shock
    P_s[i] = P_shock * (r_s[i]/r_shock)**2 * P_dimensionless(v_bar)
    v_s[i] = (4./(5.*(gamma+1.)) * r_s[i] / t * v_bar)

u_s = P_s / ((gamma - 1.)*rho_s) # internal energy
s_s = P_s / rho_s**gamma # entropic function
rho_s[0] = 0.
P_s[0] = P_s[1] # dirty...


savetxt("rho.dat", column_stack((r_s, rho_s)))
savetxt("P.dat", column_stack((r_s, P_s)))
savetxt("v.dat", column_stack((r_s, v_s)))
savetxt("u.dat", column_stack((r_s, u_s)))
savetxt("s.dat", column_stack((r_s, s_s)))



