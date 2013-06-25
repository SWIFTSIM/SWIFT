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

# Computes the analytical solution of the Gresho-Chan vortex
# The script works for a given initial box and background pressure and computes the solution for any time t (The solution is constant over time).
# The code writes five files rho.dat, P.dat, v.dat, u.dat and s.dat with the density, pressure, internal energy and
# entropic function on N radial points between r=0 and r=R_max.

# Parameters
rho0 = 1.         # Background Density
P0 = 0.           # Background Pressure
gamma = 5./3.     # Gas polytropic index
N = 1000          # Number of radial points
R_max = 1.        # Maximal radius

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

r = arange(0, R_max, R_max / N)
rho = ones(N)
P = zeros(N)
v = zeros(N)
u = zeros(N)
s = zeros(N)


for i in range(N):
    if r[i] < 0.2:
        P[i] = P0 + 5. + 12.5*r[i]**2
        v[i] = 5.*r[i]
    elif r[i] < 0.4:
        P[i] = P0 + 9. + 12.5*r[i]**2 - 20.*r[i] + 4.*log(r[i]/0.2)
        v[i] = 2. -5.*r[i]
    else:
        P[i] = P0 + 3. + 4.*log(2.)
        v[i] = 0.
    rho[i] = rho0
    s[i] = P[i] / rho[i]**gamma
    u[i] = P[i] /((gamma - 1.)*rho[i])

savetxt("rho.dat", column_stack((r, rho)))
savetxt("P.dat", column_stack((r, P)))
savetxt("v.dat", column_stack((r, v)))
savetxt("u.dat", column_stack((r, u)))
savetxt("s.dat", column_stack((r, s)))



