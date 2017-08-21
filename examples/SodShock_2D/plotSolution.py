###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Computes the analytical solution of the Sod shock and plots the SPH answer
 

# Generates the analytical  solution for the Sod shock test case
# The script works for a given left (x<0) and right (x>0) state and computes the solution at a later time t.
# This follows the solution given in (Toro, 2009)


# Parameters
gas_gamma = 5./3.      # Polytropic index
rho_L = 1.             # Density left state
rho_R = 0.140625       # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state


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
sim = h5py.File("sodShock_%04d.hdf5"%snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:,0]
v = sim["/PartType0/Velocities"][:,0]
u = sim["/PartType0/InternalEnergy"][:]
S = sim["/PartType0/Entropy"][:]
P = sim["/PartType0/Pressure"][:]
rho = sim["/PartType0/Density"][:]

N = 1000  # Number of points
x_min = -1.
x_max = 1.
x += x_min


# Bin te data
x_bin_edge = np.arange(-0.6, 0.6, 0.02)
x_bin = 0.5*(x_bin_edge[1:] + x_bin_edge[:-1])
rho_bin,_,_ = stats.binned_statistic(x, rho, statistic='mean', bins=x_bin_edge)
v_bin,_,_ = stats.binned_statistic(x, v, statistic='mean', bins=x_bin_edge)
P_bin,_,_ = stats.binned_statistic(x, P, statistic='mean', bins=x_bin_edge)
S_bin,_,_ = stats.binned_statistic(x, S, statistic='mean', bins=x_bin_edge)
u_bin,_,_ = stats.binned_statistic(x, u, statistic='mean', bins=x_bin_edge)
rho2_bin,_,_ = stats.binned_statistic(x, rho**2, statistic='mean', bins=x_bin_edge)
v2_bin,_,_ = stats.binned_statistic(x, v**2, statistic='mean', bins=x_bin_edge)
P2_bin,_,_ = stats.binned_statistic(x, P**2, statistic='mean', bins=x_bin_edge)
S2_bin,_,_ = stats.binned_statistic(x, S**2, statistic='mean', bins=x_bin_edge)
u2_bin,_,_ = stats.binned_statistic(x, u**2, statistic='mean', bins=x_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin**2)
v_sigma_bin = np.sqrt(v2_bin - v_bin**2)
P_sigma_bin = np.sqrt(P2_bin - P_bin**2)
S_sigma_bin = np.sqrt(S2_bin - S_bin**2)
u_sigma_bin = np.sqrt(u2_bin - u_bin**2)


# Analytic solution
c_L = sqrt(gas_gamma * P_L / rho_L)   # Speed of the rarefaction wave
c_R = sqrt(gas_gamma * P_R / rho_R)   # Speed of the shock front

# Helpful variable
Gama = (gas_gamma - 1.) / (gas_gamma + 1.)
beta = (gas_gamma - 1.) / (2. * gas_gamma)

# Characteristic function and its derivative, following Toro (2009)
def compute_f(P_3, P, c):
    u = P_3 / P
    if u > 1:
        term1 = gas_gamma*((gas_gamma+1.)*u + gas_gamma-1.)
        term2 = sqrt(2./term1)
        fp = (u - 1.)*c*term2
        dfdp = c*term2/P + (u - 1.)*c/term2*(-1./term1**2)*gas_gamma*(gas_gamma+1.)/P
    else:
        fp = (u**beta - 1.)*(2.*c/(gas_gamma-1.))
        dfdp = 2.*c/(gas_gamma-1.)*beta*u**(beta-1.)/P
    return (fp, dfdp)

# Solution of the Riemann problem following Toro (2009) 
def RiemannProblem(rho_L, P_L, v_L, rho_R, P_R, v_R):
    P_new = ((c_L + c_R + (v_L - v_R)*0.5*(gas_gamma-1.))/(c_L / P_L**beta + c_R / P_R**beta))**(1./beta)
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
    v_right = v_R + c_R**2*(P_3/P_R - 1.)/(gas_gamma*(v_3-v_R))
else:
    v_right = c_R + 0.5*(gas_gamma+1.)*v_3 - 0.5*(gas_gamma-1.)*v_R

if shock_L:
    v_left = v_L + c_L**2*(P_3/p_L - 1.)/(gas_gamma*(v_3-v_L))
else:
    v_left = c_L - 0.5*(gas_gamma+1.)*v_3 + 0.5*(gas_gamma-1.)*v_L

# Compute position of the transitions
x_23 = -fabs(v_left) * time
if shock_L :
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
            rho_s[i] = rho_L*(Gama + P_3/P_L)/(1. + Gama * P_3/P_L)
            P_s[i] = P_3
            v_s[i] = v_3
        else:
            rho_s[i] = rho_L*(Gama * (0. - x_s[i])/(c_L * time) + Gama * v_L/c_L + (1.-Gama))**(2./(gas_gamma-1.))
            P_s[i] = P_L*(rho_s[i] / rho_L)**gas_gamma
            v_s[i] = (1.-Gama)*(c_L -(0. - x_s[i]) / time) + Gama*v_L
    if x_s[i] >= x_23 and x_s[i] < x_34:
        if shock_L:
            rho_s[i] = rho_L*(Gama + P_3/P_L)/(1+Gama * P_3/p_L)
        else:
            rho_s[i] = rho_L*(P_3 / P_L)**(1./gas_gamma)
        P_s[i] = P_3
        v_s[i] = v_3
    if x_s[i] >= x_34 and x_s[i] < x_45:
        if shock_R:
            rho_s[i] = rho_R*(Gama + P_3/P_R)/(1. + Gama * P_3/P_R)
        else:
            rho_s[i] = rho_R*(P_3 / P_R)**(1./gas_gamma)
        P_s[i] = P_3
        v_s[i] = v_3
    if x_s[i] >= x_45 and x_s[i] < x_56:
        if shock_R:
            rho_s[i] = rho_R
            P_s[i] = P_R
            v_s[i] = v_R
        else:
            rho_s[i] = rho_R*(Gama*(x_s[i])/(c_R*time) - Gama*v_R/c_R + (1.-Gama))**(2./(gas_gamma-1.))
            P_s[i] = p_R*(rho_s[i]/rho_R)**gas_gamma
            v_s[i] = (1.-Gama)*(-c_R - (-x_s[i])/time) + Gama*v_R
    if x_s[i] >= x_56:
        rho_s[i] = rho_R
        P_s[i] = P_R
        v_s[i] = v_R


# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.))  #internal energy
s_s = P_s / rho_s**gas_gamma # entropic function
        

# Plot the interesting quantities
figure()

# Velocity profile --------------------------------
subplot(231)
plot(x, v, '.', color='r', ms=0.2)
plot(x_s, v_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(x_bin, v_bin, yerr=v_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)
xlim(-0.5, 0.5)
ylim(-0.1, 0.95)

# Density profile --------------------------------
subplot(232)
plot(x, rho, '.', color='r', ms=0.2)
plot(x_s, rho_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(x_bin, rho_bin, yerr=rho_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.05, 1.1)

# Pressure profile --------------------------------
subplot(233)
plot(x, P, '.', color='r', ms=0.2)
plot(x_s, P_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(x_bin, P_bin, yerr=P_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.01, 1.1)

# Internal energy profile -------------------------
subplot(234)
plot(x, u, '.', color='r', ms=0.2)
plot(x_s, u_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(x_bin, u_bin, yerr=u_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.8, 2.2)

# Entropy profile ---------------------------------
subplot(235)
plot(x, S, '.', color='r', ms=0.2)
plot(x_s, s_s, '--', color='k', alpha=0.8, lw=1.2)
errorbar(x_bin, S_bin, yerr=S_sigma_bin, fmt='.', ms=8.0, color='b', lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)
xlim(-0.5, 0.5)
ylim(0.8, 3.8)

# Information -------------------------------------
subplot(236, frameon=False)

text(-0.49, 0.9, "Sod shock with  $\\gamma=%.3f$ in 2D at $t=%.2f$"%(gas_gamma,time), fontsize=10)
text(-0.49, 0.8, "Left:~~ $(P_L, \\rho_L, v_L) = (%.3f, %.3f, %.3f)$"%(P_L, rho_L, v_L), fontsize=10)
text(-0.49, 0.7, "Right: $(P_R, \\rho_R, v_R) = (%.3f, %.3f, %.3f)$"%(P_R, rho_R, v_R), fontsize=10)
plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
text(-0.49, 0.5, "$\\textsc{Swift}$ %s"%git, fontsize=10)
text(-0.49, 0.4, scheme, fontsize=10)
text(-0.49, 0.3, kernel, fontsize=10)
text(-0.49, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)"%(neighbours, eta), fontsize=10)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])


savefig("SodShock.png", dpi=200)
