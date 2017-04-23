###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016  Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import integrate
import distinct_colours as colours
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from matplotlib.font_manager import FontProperties
import numpy
import math

params = {'axes.labelsize': 9,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 10,
'xtick.labelsize': 9,
'ytick.labelsize': 9,
'text.usetex': True,
'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.115,
'figure.subplot.right'   : 0.99  ,
'figure.subplot.bottom'  : 0.08  ,
'figure.subplot.top'     : 0.99  ,
'figure.subplot.wspace'  : 0.  ,
'figure.subplot.hspace'  : 0.  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})

# Parameters
epsilon = 2.
epsilon_plummer = epsilon 
r_min = 0.
r_max = 4
r_max_plot = 2.6

# Radius
r = linspace(r_min, r_max, 401)
r[0] += 1e-9
u = r / epsilon

# Newtonian solution
phi_newton = 1. / r
F_newton = 1. / r**2

# Softened potential
phi = np.zeros(np.size(r))
W = np.zeros(np.size(r))
F = np.zeros(np.size(r))
phi_plummer = np.zeros(np.size(r))
F_plummer = np.zeros(np.size(r))
for i in range(np.size(r)):
    if r[i] > epsilon:
        phi[i] = 1. / r[i]
        W[i] = 0.
        F[i] = 1. / r[i]**2
    else:
        phi[i] = (-1./epsilon) * (3.*u[i]**7 - 15.*u[i]**6 + 28.*u[i]**5 - 21.*u[i]**4 + 7.*u[i]**2 - 3.)
        W[i] = (21. / 2.*math.pi) * (4.*u[i]**5 - 15.*u[i]**4 + 20.*u[i]**3 - 10.*u[i]**2 + 1.) / epsilon**6
        F[i] = (1./epsilon**2) * (21.*u[i]**6 - 90*u[i]**5 + 140.*u[i]**4 - 84.*u[i]**3 + 14*u[i])
        
    phi_plummer[i] = (1./epsilon_plummer) * (1 + (u[i])**2)**(-1./2.)
    F_plummer[i] = (1./epsilon_plummer**3) * r[i] / (1 + (u[i])**2)**(3./2.)
        
figure()

# Potential
subplot(311)
plot(r, phi_newton, 'g--', lw=0.8, label="${\\rm Newtonian}$")
plot(r, phi_plummer, 'b-.', lw=0.8, label="${\\rm Plummer}$")
plot(r, phi, 'r-', lw=1, label="${\\rm SWIFT}$")
plot([epsilon, epsilon], [0, 10], 'k:', alpha=0.5, lw=0.5)
legend(loc="upper right", frameon=False, handletextpad=0.1, handlelength=3.2, fontsize=8)
ylim(0, 2.1)
ylabel("$\\phi(r)$", labelpad=0)


xlim(0,r_max_plot)
xticks([0., 0.5, 1., 1.5, 2., 2.5], ["", "", "", "", "", ""])

# Force
subplot(312)
plot(r, F_newton, 'g--', lw=0.8)
plot(r, F_plummer, 'b-.', lw=0.8)
plot(r, F, 'r-', lw=1)
plot([epsilon, epsilon], [0, 10], 'k:', alpha=0.5, lw=0.5)

xlim(0,r_max_plot)
xticks([0., 0.5, 1., 1.5, 2., 2.5], ["", "", "", "", "", ""])

ylim(0, 0.95)
ylabel("$|\\mathbf{\\nabla}\\phi(r)|$", labelpad=0)

# Density
subplot(313)
plot(r, W, 'r-', lw=1)
plot([epsilon, epsilon], [0, 10], 'k:', alpha=0.5, lw=0.5)
xlim(0,r_max_plot)
xticks([0., 0.5, 1., 1.5, 2., 2.5], ["$0$", "$0.5$", "$1$", "$1.5$", "$2$", "$2.5$"])
xlabel("$r$", labelpad=-1)

ylim(0., 0.58)
ylabel("$\\rho(r)$", labelpad=0)

savefig("potential.pdf")




#Construct potential
# phi = np.zeros(np.size(r))
# for i in range(np.size(r)):
#     if r[i] > 2*epsilon:
#         phi[i] = 1./ r[i]
#     elif r[i] > epsilon:
#         phi[i] = -(1./epsilon) * ((32./3.)*u[i]**2 - (48./3.)*u[i]**3 + (38.4/4.)*u[i]**4 - (32./15.)*u[i]**5 + (2./30.)*u[i]**(-1) - (9/5.))
#     else:
#         phi[i] = -(1./epsilon) * ((32./6.)*u[i]**2 - (38.4/4.)*u[i]**4 + (32./5.)*u[i]**4 - (7./5.))
