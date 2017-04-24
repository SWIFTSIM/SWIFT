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
W_newton = 0. * r

# Softened potential
phi = np.zeros(np.size(r))
W = np.zeros(np.size(r))
F = np.zeros(np.size(r))
for i in range(np.size(r)):
    if r[i] > epsilon:
        phi[i] = 1. / r[i]
        W[i] = 0.
        F[i] = 1. / r[i]**2
    else:
        phi[i] = (-1./epsilon) * (3.*u[i]**7 - 15.*u[i]**6 + 28.*u[i]**5 - 21.*u[i]**4 + 7.*u[i]**2 - 3.)
        W[i] = (21. / (2.*math.pi)) * (4.*u[i]**5 - 15.*u[i]**4 + 20.*u[i]**3 - 10.*u[i]**2 + 1.) / epsilon**3
        F[i] = (1./epsilon**2) * (21.*u[i]**6 - 90*u[i]**5 + 140.*u[i]**4 - 84.*u[i]**3 + 14*u[i])

plummer_equivalent_factor = phi[0] * epsilon

print "Plummer-equivalent factor:", plummer_equivalent_factor
        
epsilon_plummer = epsilon / plummer_equivalent_factor
        
# Plummer potential
phi_plummer = (1. / epsilon_plummer) * (1 + (r / epsilon_plummer)**2)**(-1./2.)
F_plummer = (1. / epsilon_plummer**3) * r / (1 + (r / epsilon_plummer )**2)**(3./2.)
def eta_plummer(r):
    return (3. / (4.*math.pi)) * 1. / (1 + r**2)**(5./2.)
W_plummer = (1. / epsilon_plummer**3) * eta_plummer(r / epsilon_plummer)


# Gadget-2 potential
epsilon_gadget = epsilon #/ plummer_equivalent_factor * 2.8
phi_gadget2 = np.zeros(np.size(r))
W_gadget2 = np.zeros(np.size(r))
F_gadget2 = np.zeros(np.size(r))
for i in range(np.size(r)):
    if r[i] > epsilon_gadget:
        phi_gadget2[i] = 1. / r[i]
        W_gadget2[i] = 0.
        F_gadget2[i] = 1. / r[i]**2
    elif r[i] > epsilon_gadget / 2.:
        phi_gadget2[i] = -((32/3.)*u[i]**2 - 16.*u[i]**3 + (96./10.)*u[i]**4 - (64./30.)*u[i]**5 - (16./5.) + 1./(15.*u[i])  )/ (epsilon_gadget)
        W_gadget2[i] = (8. / math.pi) * (2. * (1- u[i])**3) / epsilon_gadget**3
        F_gadget2[i] = u[i] * (21.333333 - 48*u[i] + 38.4*u[i]**2 - 10.6666667*u[i]**3 - 0.06666667*u[i]**-3) / epsilon_gadget**2   
    else:
        phi_gadget2[i] = -((16./3.)*u[i]**2 - (96./10.)*u[i]**4 + (64./10.)*u[i]**5 - (14./5.)  ) / (epsilon_gadget)
        W_gadget2[i] = (8. / math.pi) * (1. - 6.*u[i]**2 + 6.*u[i]**3) / epsilon_gadget**3
        F_gadget2[i] = u[i] * (10.666667 + u[i]**2 * (32. * u[i] - 38.4)) / epsilon_gadget**2

figure()

# Potential
subplot(311)
plot(r, phi_newton, 'g--', lw=0.8, label="${\\rm Newtonian}$")
plot(r, phi_plummer, 'b:', lw=0.8, label="${\\rm Plummer}$")
plot(r, phi_gadget2, 'm-.', lw=0.8, label="${\\rm Gadget-2}$")
plot(r, phi, 'r-', lw=0.8, label="${\\rm SWIFT}$")
plot([epsilon, epsilon], [-10, 10], 'k-', alpha=0.5, lw=0.5)
#text(epsilon+0.05, 2., "$\\epsilon$", color='k', alpha=0.5, rotation=90, va="top", ha="left", fontsize=8) 
text(1.2, 1.6, "$\\epsilon_{\\rm{Plummer}} = \\frac{\\epsilon}{%d}$"%plummer_equivalent_factor, color='k', fontsize=8, backgroundcolor='w')
ylim(0, 2.1)
ylabel("$|\\phi(r)|$", labelpad=0)


xlim(0,r_max_plot)
xticks([0., 0.5, 1., 1.5, 2., 2.5], ["", "", "", "", "", ""])

# Force
subplot(312)
plot(r, F_newton, 'g--', lw=0.8)
plot(r, F_plummer, 'b:', lw=0.8)
plot(r, F_gadget2, 'm-.', lw=0.8)
plot(r, F, 'r-', lw=0.8)
plot([epsilon, epsilon], [0, 10], 'k-', alpha=0.5, lw=0.5)

xlim(0,r_max_plot)
xticks([0., 0.5, 1., 1.5, 2., 2.5], ["", "", "", "", "", ""])

ylim(0, 0.95)
ylabel("$|\\overrightarrow{\\nabla}\\phi(r)|$", labelpad=0)

# Density
subplot(313)
plot(r, W_newton, 'g--', lw=0.8, label="${\\rm Newtonian}$")
plot(r, W_plummer, 'b:', lw=0.8, label="${\\rm Plummer}$")
plot(r, W_gadget2, 'm-.', lw=0.8, label="${\\rm Gadget\\textendash 2}$")
plot(r, W, 'r-', lw=0.8, label="${\\rm SWIFT}$")
plot([epsilon, epsilon], [0, 10], 'k-', alpha=0.5, lw=0.5)
xlim(0,r_max_plot)
xticks([0., 0.5, 1., 1.5, 2., 2.5], ["$%.1f$"%(0./epsilon), "", "$%.1f$"%(1./epsilon), "", "$%.1f$"%(2./epsilon)])
xlabel("$r/\\epsilon$", labelpad=-4)
legend(loc="upper right", frameon=True, handletextpad=0.1, handlelength=3.2, fontsize=8)

ylim(0., 0.92)
yticks([0, 0.2, 0.4, 0.6, 0.8], ["$0$", "$0.2$", "$0.4$", "$0.6$", "$0.8$"])
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
