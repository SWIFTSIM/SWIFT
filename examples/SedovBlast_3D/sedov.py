###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
import scipy.integrate as integrate
import scipy.optimize as optimize
import os

# Calculate the analytical solution of the Sedov-Taylor shock wave for a given
# number of dimensions, gamma and time. We assume dimensionless units and a
# setup with constant density 1, pressure and velocity 0. An energy 1 is
# inserted in the center at t 0.
#
# The solution is a self-similar shock wave, which was described in detail by
# Sedov (1959). We follow his notations and solution method.
#
# The position of the shock at time t is given by
#  r2 = (E/rho1)^(1/(2+nu)) * t^(2/(2+nu))
# the energy E is related to the inserted energy E0 by E0 = alpha*E, with alpha
# a constant which has to be calculated by imposing energy conservation.
#
# The density for a given radius at a certain time is determined by introducing
# the dimensionless position coordinate lambda = r/r2. The density profile as
# a function of lambda is constant in time and given by
#  rho = rho1 * R(V, gamma, nu)
# and
#  V = V(lambda, gamma, nu)
#
# The function V(lambda, gamma, nu) is found by solving a differential equation
# described in detail in Sedov (1959) chapter 4, section 5. Alpha is calculated
# from the integrals in section 11 of the same chapter.
#
# Numerically, the complete solution requires the use of 3 quadratures and 1
# root solver, which are implemented using the GNU Scientific Library (GSL).
# Since some quadratures call functions that themselves contain quadratures,
# the problem is highly non-linear and complex and takes a while to converge.
# Therefore, we tabulate the alpha values and profile values the first time
# a given set of gamma and nu values is requested and reuse these tabulated
# values.
#
# Reference:
#  Sedov (1959): Sedov, L., Similitude et dimensions en mecanique (7th ed.;
#                Moscou: Editions Mir) - french translation of the original
#                book from 1959.

# dimensionless variable z = gamma*P/R as a function of V (and gamma and nu)
# R is a dimensionless density, while P is a dimensionless pressure
# z is hence a sort of dimensionless sound speed
# The expression below corresponds to eq. 11.9 in Sedov (1959), chapter 4
def _z(V, gamma, nu):
    if V == 2./(nu+2.)/gamma:
        return 0.
    else:
        return (gamma-1.)*V*V*(V-2./(2.+nu))/2./(2./(2.+nu)/gamma-V)

# differential equation that needs to be solved to obtain lambda for a given V
# corresponds to eq. 5.11 in Sedov (1959), chapter 4 (omega = 0)
def _dlnlambda_dV(V, gamma, nu):
    nom = _z(V, gamma, nu) - (V-2./(nu+2.))*(V-2./(nu+2.))
    denom = V*(V-1.)*(V-2./(nu+2.))+nu*(2./(nu+2.)/gamma-V)*_z(V, gamma, nu)
    return nom/denom

# dimensionless variable lambda = r/r2 as a function of V (and gamma and nu)
# found by solving differential equation 5.11 in Sedov (1959), chapter 4
# (omega = 0)
def _lambda(V, gamma, nu):
    if V == 2./(nu+2.)/gamma:
        return 0.
    else:
        V0 = 4./(nu+2.)/(gamma+1.)
        integral, err = integrate.quad(_dlnlambda_dV, V, V0, (gamma, nu),
                                       limit = 8000)
        return np.exp(-integral)

# dimensionless variable R = rho/rho1 as a function of V (and gamma and nu)
# found by inverting eq. 5.12 in Sedov (1959), chapter 4 (omega = 0)
# the integration constant C is found by inserting the R, V and z values
# at the shock wave, where lambda is 1. These correspond to eq. 11.8 in Sedov
# (1959), chapter 4.
def _R(V, gamma, nu):
    if V == 2./(nu+2.)/gamma:
        return 0.
    else:
        C = 8.*gamma*(gamma-1.)/(nu+2.)/(nu+2.)/(gamma+1.)/(gamma+1.) \
                *((gamma-1.)/(gamma+1.))**(gamma-2.) \
                *(4./(nu+2.)/(gamma+1.)-2./(nu+2.))
        lambda1 = _lambda(V, gamma, nu)
        lambda5 = lambda1**(nu+2)
        return (_z(V, gamma, nu)*(V-2./(nu+2.))*lambda5/C)**(1./(gamma-2.))

# function of which we need to find the zero point to invert lambda(V)
def _lambda_min_lambda(V, lambdax, gamma, nu):
    return _lambda(V, gamma, nu) - lambdax

# dimensionless variable V = v*t/r as a function of lambda (and gamma and nu)
# found by inverting the function lamdba(V) which is found by solving
# differential equation 5.11 in Sedov (1959), chapter 4 (omega = 0)
# the invertion is done by searching the zero point of the function
# lambda_min_lambda defined above
def _V_inv(lambdax, gamma, nu):
    if lambdax == 0.:
        return 2./(2.+nu)/gamma;
    else:
        return optimize.brentq(_lambda_min_lambda, 2./(nu+2.)/gamma, 
                               4./(nu+2.)/(gamma+1.), (lambdax, gamma, nu))

# integrand of the first integral in eq. 11.24 in Sedov (1959), chapter 4
def _integrandum1(lambdax, gamma, nu):
    V = _V_inv(lambdax, gamma, nu)
    if nu == 2:
        return _R(V, gamma, nu)*V**2*lambdax**3
    else:
        return _R(V, gamma, nu)*V**2*lambdax**4

# integrand of the second integral in eq. 11.24 in Sedov (1959), chapter 4
def _integrandum2(lambdax, gamma, nu):
    V = _V_inv(lambdax, gamma, nu)
    if V == 2./(nu+2.)/gamma:
        P = 0.
    else:
        P = _z(V, gamma, nu)*_R(V, gamma, nu)/gamma
    if nu == 2:
        return P*lambdax**3
    else:
        return P*lambdax**4

# calculate alpha = E0/E
# this corresponds to eq. 11.24 in Sedov (1959), chapter 4
def get_alpha(gamma, nu):
  integral1, err1 = integrate.quad(_integrandum1, 0., 1., (gamma, nu))
  integral2, err2 = integrate.quad(_integrandum2, 0., 1., (gamma, nu))
  
  if nu == 2:
    return np.pi*integral1+2.*np.pi/(gamma-1.)*integral2
  else:
    return 2.*np.pi*integral1+4.*np.pi/(gamma-1.)*integral2

# get the analytical solution for the Sedov-Taylor blastwave given an input
# energy E, adiabatic index gamma, and number of dimensions nu, at time t and
# with a maximal outer region radius maxr
def get_analytical_solution(E, gamma, nu, t, maxr = 1.):
    # we check for the existance of a datafile with precomputed alpha and
    # profile values
    # if it does not exist, we calculate it here and write it out
    # calculation of alpha and the profile takes a while...
    lvec = np.zeros(1000)
    Rvec = np.zeros(1000)
    fname = "sedov_profile_gamma_{gamma}_nu_{nu}.dat".format(gamma = gamma,
                                                             nu = nu)
    if os.path.exists(fname):
        file = open(fname, "r")
        lines = file.readlines()
        alpha = float(lines[0])
        for i in range(len(lines)-1):
            data = lines[i+1].split()
            lvec[i] = float(data[0])
            Rvec[i] = float(data[1])
    else:
        alpha = get_alpha(gamma, nu)
        for i in range(1000):
            lvec[i] = (i+1)*0.001
            V = _V_inv(lvec[i], gamma, nu)
            Rvec[i] = _R(V, gamma, nu)
        file = open(fname, "w")
        file.write("{alpha}\n".format(alpha = alpha))
        for i in range(1000):
            file.write("{l}\t{R}\n".format(l = lvec[i], R = Rvec[i]))

    xvec = np.zeros(1002)
    rhovec = np.zeros(1002)
    if nu == 2:
        r2 = (E/alpha)**0.25*np.sqrt(t)
    else:
        r2 = (E/alpha)**0.2*t**0.4

    for i in range(1000):
        xvec[i] = lvec[i]*r2
        rhovec[i] = Rvec[i]
    xvec[1000] = 1.001*r2
    xvec[1001] = maxr
    rhovec[1000] = 1.
    rhovec[1001] = 1.
    
    return xvec, rhovec

def main():
    E = 1.
    gamma = 1.66667
    nu = 3
    t = 0.001
    x, rho = get_analytical_solution(E, gamma, nu, t)
    for i in range(len(x)):
        print x[i], rho[i]

if __name__ == "__main__":
    main()
