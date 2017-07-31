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
from scipy import special
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from matplotlib.font_manager import FontProperties
import numpy
import math

params = {'axes.labelsize': 9,
'axes.titlesize': 10,
'font.size': 10,
'legend.fontsize': 10,
'xtick.labelsize': 8,
'ytick.labelsize': 8,
'text.usetex': True,
'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.12,
'figure.subplot.right'   : 0.99  ,
'figure.subplot.bottom'  : 0.09  ,
'figure.subplot.top'     : 0.99  ,
'figure.subplot.wspace'  : 0.  ,
'figure.subplot.hspace'  : 0.  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})
colors=['#4477AA', '#CC6677', '#DDCC77', '#117733']


# Parameters
r_s = 2.
r_min = 1e-2
r_max = 1.5e2

# Radius
r = logspace(log10(r_min), log10(r_max), 401)
r_rs = r / r_s

k = logspace(log10(r_min/r_s**2), log10(r_max/r_s**2), 401)
k_rs = k * r_s

# Newtonian solution
phi_newton = 1. / r
phit_newton = 1. / k**2
force_newton = 1. / r**2

def my_exp(x):
    return 1. + x + (x**2 / 2.) + (x**3 / 6.) + (x**4 / 24.) + (x**5 / 120.)# + (x**6 / 720.)
    #return exp(x)

def csch(x): # hyperbolic cosecant
    return 1. / sinh(x)

def sigmoid(x):
    return my_exp(x) / (my_exp(x) + 1.)

def d_sigmoid(x):
    return my_exp(x) / ((my_exp(x) + 1)**2)

def swift_corr(x):
    return 2 * sigmoid( 4 * x ) - 1

#figure()
#x = linspace(-4, 4, 100)
#plot(x, special.erf(x), '-', color=colors[0])
#plot(x, swift_corr(x), '-', color=colors[1])
#plot(x, x, '-', color=colors[2])
#ylim(-1.1, 1.1)
#xlim(-4.1, 4.1)
#savefig("temp.pdf")

# Correction in real space
corr_short_gadget2 = special.erf(r / (2.*r_s))
corr_short_swift = swift_corr(r / (2.*r_s)) 
eta_short_gadget2 = special.erfc(r / 2.*r_s) + (r / (r_s * math.sqrt(math.pi))) * exp(-r**2 / (4.*r_s**2))
eta_short_swift = 4. * (r / r_s) * d_sigmoid(2. * r / r_s) - 2. * sigmoid(2 * r / r_s) + 2.

# Corection in Fourier space
corr_long_gadget2 = exp(-k**2*r_s**2)
corr_long_swift = math.pi * k * r_s * csch(0.5 * math.pi * r_s * k) / 2.

# Shortrange term
phi_short_gadget2 = (1.  / r ) * (1. - corr_short_gadget2)
phi_short_swift = (1.  / r ) * (1. - corr_short_swift)
force_short_gadget2 = (1. / r**2) * eta_short_gadget2
force_short_swift = (1. / r**2) * eta_short_swift

# Long-range term
phi_long_gadget2 = (1.  / r ) * corr_short_gadget2
phi_long_swift = (1.  / r ) * corr_short_swift
phit_long_gadget2 = corr_long_gadget2 / k**2
phit_long_swift = corr_long_swift / k**2




figure()

# Potential
subplot(311, xscale="log", yscale="log")

plot(r_rs, phi_newton, '--', lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(r_rs, phi_short_gadget2, '-', lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(r_rs, phi_short_swift, '-', lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)

xlim(1.1*r_min/ r_s, 0.9*r_max / r_s)
ylim(1.1/r_max, 0.9/r_min)
ylabel("$\\varphi_s(r)$", labelpad=-3)

legend(loc="upper right", frameon=True, handletextpad=0.1, handlelength=3.2, fontsize=8)

# Correction
subplot(312, xscale="log", yscale="log")
plot(r_rs, np.ones(np.size(r)), '--', lw=1.4, color=colors[0])
plot(r_rs, 1. - corr_short_gadget2, '-', lw=1.4, color=colors[2])
plot(r_rs, 1. - corr_short_swift, '-', lw=1.4, color=colors[3])
plot(r_rs, np.ones(np.size(r))*0.01, 'k:', alpha=0.5, lw=0.5)
plot([1., 1.], [-1e5, 1e5], 'k-', alpha=0.5, lw=0.5)

yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlim(1.1*r_min/r_s, 0.9*r_max/r_s)
ylim(3e-3, 1.5)
#ylabel("$\\chi_s(r)$", labelpad=-3)
ylabel("$\\varphi_s(r) \\times r$", labelpad=-2)

# 1 - Correction
subplot(313, xscale="log", yscale="log")
plot(r_rs, corr_short_gadget2, '-', lw=1.4, color=colors[2])
plot(r_rs, corr_short_swift, '-', lw=1.4, color=colors[3])

plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)), 'k:', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r))*0.01, 'k:', alpha=0.5, lw=0.5)

xlim(1.1*r_min/r_s, 0.9*r_max/r_s)
ylim(3e-3, 1.5)
#ylabel("$1 - \\chi_s(r)$", labelpad=-2)
ylabel("$1 - \\varphi_s(r) \\times r$", labelpad=-2)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlabel("$r / r_s$", labelpad=-3)

savefig("potential_short.pdf")

##################################################################################################


# Force
figure()
subplot(311, xscale="log", yscale="log")

plot(r_rs, force_newton, '--', lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(r_rs, force_short_gadget2, '-', lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(r_rs, force_short_swift, '-', lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)

xlim(1.1*r_min/ r_s, 0.9*r_max / r_s)
ylim(1.1/r_max**2, 0.9/r_min**2)
ylabel("$|\\mathbf{f}_s(r)|$", labelpad=-3)
yticks([1e-4, 1e-2, 1e0, 1e2], ["$10^{-4}$", "$10^{-2}$", "$10^{0}$", "$10^{2}$"])

legend(loc="upper right", frameon=True, handletextpad=0.1, handlelength=3.2, fontsize=8)

# Correction
subplot(312, xscale="log", yscale="log")
plot(r_rs, np.ones(np.size(r)), '--', lw=1.4, color=colors[0])
plot(r_rs, eta_short_gadget2, '-', lw=1.4, color=colors[2])
plot(r_rs, eta_short_swift, '-', lw=1.4, color=colors[3])
plot(r_rs, np.ones(np.size(r))*0.01, 'k:', alpha=0.5, lw=0.5)
plot([1., 1.], [-1e5, 1e5], 'k-', alpha=0.5, lw=0.5)

yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlim(1.1*r_min/r_s, 0.9*r_max/r_s)
ylim(3e-3, 1.5)
#ylabel("$\\eta_s(r)$", labelpad=-3)
ylabel("$|\\mathbf{f}_s(r)|\\times r^2$", labelpad=-2)

# 1 - Correction
subplot(313, xscale="log", yscale="log")
plot(r_rs, 1. - eta_short_gadget2, '-', lw=1.4, color=colors[2])
plot(r_rs, 1. - eta_short_swift, '-', lw=1.4, color=colors[3])

plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)), 'k:', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r))*0.01, 'k:', alpha=0.5, lw=0.5)

xlim(1.1*r_min/r_s, 0.9*r_max/r_s)
ylim(3e-3, 1.5)
#ylabel("$1 - \\eta_s(r)$", labelpad=-2)
ylabel("$1 - |\\mathbf{f}_s(r)|\\times r^2$", labelpad=-3)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlabel("$r / r_s$", labelpad=-3)

savefig("force_short.pdf")

##################################################################################################

figure()
subplot(311, xscale="log", yscale="log")

# Potential
plot(k_rs, phit_newton, '--', lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(k_rs, phit_long_gadget2, '-', lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(k_rs, phit_long_swift, '-', lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot(k_rs, -phit_long_swift, ':', lw=1.4, color=colors[3])
plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)

legend(loc="lower left", frameon=True, handletextpad=0.1, handlelength=3.2, fontsize=8)

xlim(1.1*r_min/ r_s, 0.9*r_max / r_s)
ylim(1.1/r_max**2, 0.9/r_min**2)
ylabel("$\\tilde{\\varphi_l}(k)$", labelpad=-3)
yticks([1e-4, 1e-2, 1e0, 1e2], ["$10^{-4}$", "$10^{-2}$", "$10^{0}$", "$10^{2}$"])

subplot(312, xscale="log", yscale="log")

# Potential normalized
plot(k_rs, phit_newton * k**2, '--', lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(k_rs, phit_long_gadget2 * k**2, '-', lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(k_rs, phit_long_swift * k**2, '-', lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r))*0.01, 'k:', alpha=0.5, lw=0.5)

xlim(1.1*r_min/ r_s, 0.9*r_max / r_s)
ylim(3e-3, 1.5)
ylabel("$k^2 \\times \\tilde{\\varphi_l}(k)$", labelpad=-3)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])

subplot(313, xscale="log", yscale="log")

plot(k_rs, 1. - phit_long_gadget2 * k**2, '-', lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(k_rs, 1. - phit_long_swift * k**2, '-', lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot([1., 1.], [1e-5, 1e5], 'k-', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)), 'k:', alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r))*0.01, 'k:', alpha=0.5, lw=0.5)

xlim(1.1*r_min/ r_s, 0.9*r_max / r_s)
ylim(3e-3, 1.5)
ylabel("$1 - k^2 \\times \\tilde{\\varphi_l}(k)$", labelpad=-3)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])

xlabel("$k \\times r_s$", labelpad=0)

savefig("potential_long.pdf")
