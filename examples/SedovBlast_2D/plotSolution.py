###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

# Computes the analytical solution of the 2D Sedov blast wave.
# The script works for a given initial box and dumped energy and computes the solution at a later time t.

# Parameters
rho_0 = 1.          # Background Density
P_0 = 1.e-6         # Background Pressure
E_0 = 1.            # Energy of the explosion
gas_gamma = 5./3.   # Gas polytropic index


# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (9.90,6.45),
'figure.subplot.left'    : 0.045,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.05,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


snap = int(sys.argv[1])


# Read the simulation data
sim = h5py.File("sedov_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

pos = sim["/PartType0/Coordinates"][:,:]
x = pos[:,0] - boxSize / 2
y = pos[:,1] - boxSize / 2
vel = sim["/PartType0/Velocities"][:,:]
r = sqrt(x**2 + y**2)
v_r = (x * vel[:,0] + y * vel[:,1]) / r
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

# Bin te data
r_bin_edge = np.arange(0., 0.5, 0.01)
r_bin = 0.5*(r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin,_,_ = stats.binned_statistic(r, rho, statistic='mean', bins=r_bin_edge)
v_bin,_,_ = stats.binned_statistic(r, v_r, statistic='mean', bins=r_bin_edge)
P_bin,_,_ = stats.binned_statistic(r, P, statistic='mean', bins=r_bin_edge)
S_bin,_,_ = stats.binned_statistic(r, S, statistic='mean', bins=r_bin_edge)
u_bin,_,_ = stats.binned_statistic(r, u, statistic='mean', bins=r_bin_edge)
rho2_bin,_,_ = stats.binned_statistic(r, rho**2, statistic='mean', bins=r_bin_edge)
v2_bin,_,_ = stats.binned_statistic(r, v_r**2, statistic='mean', bins=r_bin_edge)
P2_bin,_,_ = stats.binned_statistic(r, P**2, statistic='mean', bins=r_bin_edge)
S2_bin,_,_ = stats.binned_statistic(r, S**2, statistic='mean', bins=r_bin_edge)
u2_bin,_,_ = stats.binned_statistic(r, u**2, statistic='mean', bins=r_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)

# Now, work our the solution....

from scipy.special import gamma as Gamma
from numpy import *

def calc_a(g,nu=3):
    """ 
    exponents of the polynomials of the sedov solution
    g - the polytropic gamma
    nu - the dimension
    """
    a = [0]*8
   
    a[0] = 2.0 / (nu + 2)
    a[2] = (1-g) / (2*(g-1) + nu)
    a[3] = nu / (2*(g-1) + nu)
    a[5] = 2 / (g-2)
    a[6] = g / (2*(g-1) + nu)
   
    a[1] = (((nu+2)*g)/(2.0+nu*(g-1.0)) ) * ( (2.0*nu*(2.0-g))/(g*(nu+2.0)**2) - a[2])
    a[4] = a[1]*(nu+2) / (2-g)
    a[7] = (2 + nu*(g-1))*a[1]/(nu*(2-g))
    return a

def calc_beta(v, g, nu=3):
    """ 
    beta values for the sedov solution (coefficients of the polynomials of the similarity variables) 
    v - the similarity variable
    g - the polytropic gamma
    nu- the dimension
    """

    beta = (nu+2) * (g+1) * array((0.25, (g/(g-1))*0.5,
            -(2 + nu*(g-1))/2.0 / ((nu+2)*(g+1) -2*(2 + nu*(g-1))),
     -0.5/(g-1)), dtype=float64)

    beta = outer(beta, v)

    beta += (g+1) * array((0.0,  -1.0/(g-1),
                           (nu+2) / ((nu+2)*(g+1) -2.0*(2 + nu*(g-1))),
                           1.0/(g-1)), dtype=float64).reshape((4,1))

    return beta


def sedov(t, E0, rho0, g, n=1000, nu=3):
    """ 
    solve the sedov problem
    t - the time
    E0 - the initial energy
    rho0 - the initial density
    n - number of points (10000)
    nu - the dimension
    g - the polytropic gas gamma
    """
    # the similarity variable
    v_min = 2.0 / ((nu + 2) * g)
    v_max = 4.0 / ((nu + 2) * (g + 1))

    v = v_min + arange(n) * (v_max - v_min) / (n - 1.0)

    a = calc_a(g, nu)
    beta = calc_beta(v, g=g, nu=nu)
    lbeta = log(beta)
    
    r = exp(-a[0] * lbeta[0] - a[2] * lbeta[1] - a[1] * lbeta[2])
    rho = ((g + 1.0) / (g - 1.0)) * exp(a[3] * lbeta[1] + a[5] * lbeta[3] + a[4] * lbeta[2])
    p = exp(nu * a[0] * lbeta[0] + (a[5] + 1) * lbeta[3] + (a[4] - 2 * a[1]) * lbeta[2])
    u = beta[0] * r * 4.0 / ((g + 1) * (nu + 2))
    p *= 8.0 / ((g + 1) * (nu + 2) * (nu + 2))

    # we have to take extra care at v=v_min, since this can be a special point.
    # It is not a singularity, however, the gradients of our variables (wrt v) are.
    # r -> 0, u -> 0, rho -> 0, p-> constant

    u[0] = 0.0; rho[0] = 0.0; r[0] = 0.0; p[0] = p[1]

    # volume of an n-sphere
    vol = (pi ** (nu / 2.0) / Gamma(nu / 2.0 + 1)) * power(r, nu)

    # note we choose to evaluate the integral in this way because the
    # volumes of the first few elements (i.e near v=vmin) are shrinking 
    # very slowly, so we dramatically improve the error convergence by 
    # finding the volumes exactly. This is most important for the
    # pressure integral, as this is on the order of the volume.

    # (dimensionless) energy of the model solution
    de = rho * u * u * 0.5 + p / (g - 1)
    # integrate (trapezium rule)
    q = inner(de[1:] + de[:-1], diff(vol)) * 0.5

    # the factor to convert to this particular problem
    fac = (q * (t ** nu) * rho0 / E0) ** (-1.0 / (nu + 2))

    # shock speed
    shock_speed = fac * (2.0 / (nu + 2))
    rho_s = ((g + 1) / (g - 1)) * rho0                                                                            
    r_s = shock_speed * t * (nu + 2) / 2.0
    p_s = (2.0 * rho0 * shock_speed * shock_speed) / (g + 1)
    u_s = (2.0 * shock_speed) / (g + 1)
    
    r *= fac * t
    u *= fac
    p *= fac * fac * rho0
    rho *= rho0
    return r, p, rho, u, r_s, p_s, rho_s, u_s, shock_speed


# The main properties of the solution
r_s, P_s, rho_s, v_s, r_shock, _, _, _, _ = sedov(time, E_0, rho_0, gas_gamma, 1000, 2)

# Append points for after the shock
r_s = np.insert(r_s, np.size(r_s), [r_shock, r_shock*1.5])
rho_s = np.insert(rho_s, np.size(rho_s), [rho_0, rho_0])
P_s = np.insert(P_s, np.size(P_s), [P_0, P_0])
v_s = np.insert(v_s, np.size(v_s), [0, 0])

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.))  #internal energy
s_s = P_s / rho_s**gas_gamma # entropic function



# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
plot(r, v_r, '.', color='r', ms=0.5, alpha=0.2)
plot(r_s, v_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, v_bin, yerr=v_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Radial~velocity}}~v_r$", labelpad=0)
xlim(0, 1.3 * r_shock)
ylim(-0.2, 3.8)

# Density profile --------------------------------
subplot(232)
plot(r, rho, '.', color='r', ms=0.5, alpha=0.2)
plot(r_s, rho_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=2)
xlim(0, 1.3 * r_shock)
ylim(-0.2, 5.2)

# Pressure profile --------------------------------
subplot(233)
plot(r, P, '.', color='r', ms=0.5, alpha=0.2)
plot(r_s, P_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, P_bin, yerr=P_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)
xlim(0, 1.3 * r_shock)
ylim(-1, 12.5)

# Internal energy profile -------------------------
subplot(234)
plot(r, u, '.', color='r', ms=0.5, alpha=0.2)
plot(r_s, u_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, u_bin, yerr=u_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
xlim(0, 1.3 * r_shock)
ylim(-2, 22)

# Entropy profile ---------------------------------
subplot(235)
plot(r, S, '.', color='r', ms=0.5, alpha=0.2)
plot(r_s, s_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(r_bin, S_bin, yerr=S_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Radius}}~r$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)
xlim(0, 1.3 * r_shock)
ylim(-5, 50)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Sedov blast with  $\\gamma=%.3f$ in 2D at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "Background $\\rho_0=%.2f$"%(rho_0), fontsize=10)
text(-0.49, 0.7, "Energy injected $E_0=%.2f$"%(E_0), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])


savefig("Sedov.png", dpi=200)




