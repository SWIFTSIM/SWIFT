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

params = {
    'axes.labelsize': 10,
    'axes.titlesize': 8,
    'font.size': 10,
    'legend.fontsize': 9,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.major.pad': 2.5,
    'ytick.major.pad': 2.5,
    'text.usetex': True,
'figure.figsize' : (4.15,4.15),
'figure.subplot.left'    : 0.14,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.06,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.  ,
'figure.subplot.hspace'  : 0.  ,
'lines.markersize' : 6,
'lines.linewidth' : 1.5,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


#Parameters
eta = 1.2348422195325 # Resolution (Gives 48 neighbours for a cubic spline kernel)
dx  = 1.5#4 #2.7162  # Mean inter-particle separation

#Constants
PI = math.pi

# Compute expected moothing length
h = eta * dx

# Get kernel support (Dehnen & Aly 2012, table 1) for 3D kernels
H_cubic   = 1.825742 * h
H_quartic = 2.018932 * h
H_quintic = 2.195775 * h
H_WendlandC2 = 1.936492 * h
H_WendlandC4 = 2.207940 * h
H_WendlandC6 = 2.449490 * h

# Get the number of neighbours within kernel support:
N_H_cubic = 4./3. * PI * H_cubic**3 / (dx)**3
N_H_quartic = 4./3. * PI * H_quartic**3 / (dx)**3
N_H_quintic = 4./3. * PI * H_quintic**3 / (dx)**3
N_H_WendlandC2 = 4./3. * PI * H_WendlandC2**3 / (dx)**3
N_H_WendlandC4 = 4./3. * PI * H_WendlandC4**3 / (dx)**3
N_H_WendlandC6 = 4./3. * PI * H_WendlandC6**3 / (dx)**3


print "Smoothing length: h =", h, "Cubic spline kernel support size:   H =", H_cubic, "Number of neighbours N_H =", N_H_cubic
print "Smoothing length: h =", h, "Quartic spline kernel support size: H =", H_quartic, "Number of neighbours N_H =", N_H_quartic
print "Smoothing length: h =", h, "Quintic spline kernel support size: H =", H_quintic, "Number of neighbours N_H =", N_H_quintic
print "Smoothing length: h =", h, "Wendland C2 kernel support size:    H =", H_WendlandC2, "Number of neighbours N_H =", N_H_WendlandC2
print "Smoothing length: h =", h, "Wendland C4 kernel support size:    H =", H_WendlandC4, "Number of neighbours N_H =", N_H_WendlandC4
print "Smoothing length: h =", h, "Wendland C6 kernel support size:    H =", H_WendlandC6, "Number of neighbours N_H =", N_H_WendlandC6

# Get kernel constants (Dehen & Aly 2012, table 1) for 3D kernel
C_cubic   = 16. / PI
C_quartic = 5**6 / (512 * PI)
C_quintic = 3**7 / (40 * PI)
C_WendlandC2 = 21. / (2 * PI)
C_WendlandC4 = 495. / (32 * PI)
C_WendlandC6 = 1365. / (64 * PI)

# Get the reduced kernel definitions (Dehen & Aly 2012, table 1) for 3D kernel
#def plus(u) : return maximum(0., u)
def cubic_spline(r):   return where(r > 1., 0., where(r < 0.5,
                                                      3.*r**3 - 3.*r**2 + 0.5,
                                                      -r**3 + 3.*r**2 - 3.*r + 1.) )

#return plus(1. - r)**3 - 4.*plus(1./2. - r)**3
def quartic_spline(r): return where(r > 1., 0., where(r < 0.2,
                                                      6.*r**4 - 2.4*r**2 + 46./125.,
                                                      where(r < 0.6,
                                                            -4.*r**4 + 8.*r**3  - (24./5.)*r**2 + (8./25.)*r + 44./125.,
                                                            1.*r**4 - 4.*r**3 + 6.*r**2 - 4.*r + 1.)))

#return plus(1. - r)**4 - 5.*plus(3./5. - r)**4 + 10.*plus(1./5. - r)**4
def quintic_spline(r): return where(r > 1., 0., where(r < 1./3.,
                                                      -10.*r**5 + 10.*r**4 - (20./9.)*r**2 + (22./81.),
                                                      where(r < 2./3.,
                                                            5.*r**5 - 15.*r**4 + (50./3.)*r**3 - (70./9.)*r**2 + (25./27.)*r + (17./81.),
                                                            -1.*r**5 + 5.*r**4 - 10.*r**3 + 10.*r**2 - 5.*r + 1.)))
                                                            
#return plus(1. - r)**5 - 6.*plus(2./3. - r)**5 + 15.*plus(1./3. - r)**5
def wendlandC2(r):     return where(r > 1., 0., 4.*r**5 - 15.*r**4 + 20.*r**3 - 10*r**2 + 1.)
def wendlandC4(r):     return where(r > 1., 0.,  (35./3.)*r**8 - 64.*r**7 + 140.*r**6 - (448./3.)*r**5 + 70.*r**4 - (28. /3.)*r**2 + 1.)
def wendlandC6(r):     return where(r > 1., 0., 32.*r**11 - 231.*r**10 + 704.*r**9 - 1155.*r**8 + 1056.*r**7 - 462.*r**6 + 66.*r**4 - 11.*r**2 + 1.)
def Gaussian(r,h): return (1./(0.5*pi*h**2)**(3./2.)) * exp(- 2.*r**2 / (h**2))


# Kernel definitions (3D)
def W_cubic_spline(r):   return C_cubic      * cubic_spline(r / H_cubic)     / H_cubic**3
def W_quartic_spline(r): return C_quartic    * quartic_spline(r / H_quartic) / H_quartic**3
def W_quintic_spline(r): return C_quintic    * quintic_spline(r / H_quintic) / H_quintic**3
def W_WendlandC2(r):     return C_WendlandC2 * wendlandC2(r / H_WendlandC2)  / H_WendlandC2**3
def W_WendlandC4(r):     return C_WendlandC4 * wendlandC4(r / H_WendlandC4)  / H_WendlandC4**3
def W_WendlandC6(r):     return C_WendlandC6 * wendlandC6(r / H_WendlandC6)  / H_WendlandC6**3

# PLOT STUFF
figure()
subplot(211)
xx = linspace(0., 5*h, 1000)
maxY = 1.2*Gaussian(0, h)

# Plot the kernels
plot(xx, Gaussian(xx, h), 'k-', linewidth=0.7, label="${\\rm %14s\\quad H=\\infty}$"%("Gaussian~~~~~~"))
plot(xx, W_cubic_spline(xx), 'b-', label="${\\rm %14s\\quad H=%4.3f}$"%("Cubic~spline~~", H_cubic))
plot(xx, W_quartic_spline(xx), 'c-', label="${\\rm %14s\\quad H=%4.3f}$"%("Quartic~spline", H_quartic))
plot(xx, W_quintic_spline(xx), 'g-', label="${\\rm %14s\\quad H=%4.3f}$"%("Quintic~spline", H_quintic))
plot(xx, W_WendlandC2(xx), 'r-', label="${\\rm %14s\\quad H=%4.3f}$"%("Wendland~C2~", H_WendlandC2))
plot(xx, W_WendlandC4(xx), 'm-', label="${\\rm %14s\\quad H=%4.3f}$"%("Wendland~C4~", H_WendlandC4))
plot(xx, W_WendlandC6(xx), 'y-', label="${\\rm %14s\\quad H=%4.3f}$"%("Wendland~C6~", H_WendlandC6))

# Indicate the position of H
arrow(H_cubic, 0.12*maxY , 0., -0.12*maxY*0.9, fc='b', ec='b', length_includes_head=True, head_width=0.03, head_length=0.12*maxY*0.3)
arrow(H_quartic, 0.12*maxY , 0., -0.12*maxY*0.9, fc='c', ec='c', length_includes_head=True, head_width=0.03, head_length=0.12*maxY*0.3)
arrow(H_quintic, 0.12*maxY , 0., -0.12*maxY*0.9, fc='g', ec='g', length_includes_head=True, head_width=0.03, head_length=0.12*maxY*0.3)
arrow(H_WendlandC2, 0.12*maxY , 0., -0.12*maxY*0.9, fc='r', ec='r', length_includes_head=True, head_width=0.03, head_length=0.12*maxY*0.3)
arrow(H_WendlandC4, 0.12*maxY , 0., -0.12*maxY*0.9, fc='m', ec='m', length_includes_head=True, head_width=0.03, head_length=0.12*maxY*0.3)
arrow(H_WendlandC6, 0.12*maxY , 0., -0.12*maxY*0.9, fc='y', ec='y', length_includes_head=True, head_width=0.03, head_length=0.12*maxY*0.3)

# Show h
plot([h, h], [0., maxY], 'k:', linewidth=0.5)
text(h, maxY*0.35, "$h\\equiv\\eta\\langle x\\rangle = %.4f$"%h, rotation=90, backgroundcolor='w', ha='center', va='bottom')

# Show sigma
plot([h/2, h/2], [0., maxY], 'k:', linewidth=0.5)
text(h/2, maxY*0.05, "$\\sigma\\equiv h/2$", rotation=90, backgroundcolor='w', ha='center', va='bottom')

# Show <x>
plot([dx, dx], [0., maxY], 'k:', linewidth=0.5)
text(dx, maxY*0.35, "$\\langle x\\rangle = %.1f$"%dx, rotation=90, backgroundcolor='w', ha='center', va='bottom')

xlim(0., 2.5*h)
ylim(0., maxY)
gca().xaxis.set_ticklabels([])
ylabel("$W(r,h)$", labelpad=1.5)
legend(loc="upper right", handlelength=1.2, handletextpad=0.2)


# Same but now in log space
subplot(212, yscale="log")
plot(xx, Gaussian(xx, h), 'k-', linewidth=0.7, label="${\\rm Gaussian}$")
plot(xx, W_cubic_spline(xx), 'b-', label="${\\rm Cubic~spline}$")
plot(xx, W_quartic_spline(xx), 'c-', label="${\\rm Quartic~spline}$")
plot(xx, W_quintic_spline(xx), 'g-', label="${\\rm Quintic~spline}$")
plot(xx, W_WendlandC2(xx), 'r-', label="${\\rm Wendland~C2}$")
plot(xx, W_WendlandC4(xx), 'm-', label="${\\rm Wendland~C4}$")
plot(xx, W_WendlandC6(xx), 'y-', label="${\\rm Wendland~C6}$")

# Show h
plot([h, h], [0., 1.], 'k:', linewidth=0.5)

# Show sigma
plot([h/2, h/2], [0., 1.], 'k:', linewidth=0.5)

# Show <x>
plot([dx, dx], [0., 1.], 'k:', linewidth=0.5)


# Show plot properties
text(h/5., 1e-3, "$\\langle x \\rangle = %3.1f$"%(dx), va="top", backgroundcolor='w')
text(h/5.+0.06, 3e-4, "$\\eta = %5.4f$"%(eta), va="top", backgroundcolor='w')

# Show number of neighbours
text(1.9*h, 2e-1/2.9**0, "$N_{\\rm ngb}=\\infty$", fontsize=10)
text(1.9*h, 2e-1/2.9**1, "$N_{\\rm ngb}=%3.1f$"%(N_H_cubic), color='b', fontsize=9)
text(1.9*h, 2e-1/2.9**2, "$N_{\\rm ngb}=%3.1f$"%(N_H_quartic), color='c', fontsize=9)
text(1.9*h, 2e-1/2.9**3, "$N_{\\rm ngb}=%3.1f$"%(N_H_quintic), color='g', fontsize=9)
text(1.9*h, 2e-1/2.9**4, "$N_{\\rm ngb}=%3.1f$"%(N_H_WendlandC2), color='r', fontsize=9)
text(1.9*h, 2e-1/2.9**5, "$N_{\\rm ngb}=%3.1f$"%(N_H_WendlandC4), color='m', fontsize=9)
text(1.9*h, 2e-1/2.9**6, "$N_{\\rm ngb}=%3.0f$"%(N_H_WendlandC6), color='y', fontsize=9)

xlim(0., 2.5*h)
ylim(1e-5, 0.7)
xlabel("$r$", labelpad=0)
ylabel("$W(r,h)$", labelpad=0.5)

savefig("kernels.pdf")




################################
# Now, let's work on derivatives
################################

# Get the derivative of the reduced kernel definitions for 3D kernels
def d_cubic_spline(r):   return where(r > 1., 0., where(r < 0.5,
                                                        9.*r**2 - 6.*r,
                                                        -3.*r**2 + 6.*r - 3.) )

def d_quartic_spline(r): return where(r > 1., 0., where(r < 0.2,
                                                        24.*r**3 - 4.8*r,
                                                        where(r < 0.6,
                                                              -16.*r**3 + 24.*r**2  - (48./5.)*r + (8./25.),
                                                              4.*r**3 - 12.*r**2 + 12.*r - 4.)))

def d_quintic_spline(r): return where(r > 1., 0., where(r < 1./3.,
                                                        -50.*r**4 + 40.*r**3 - (40./9.)*r,
                                                        where(r < 2./3.,
                                                              25.*r**4 - 60.*r**3 + 50.*r**2 - (140./9.)*r + (25./27.),
                                                              -5.*r**4 + 20.*r**3 - 30.*r**2 + 20.*r - 5.)))

def d_wendlandC2(r):     return where(r > 1., 0., 20.*r**4 - 60.*r**3 + 60.*r**2 - 20.*r)
def d_wendlandC4(r):     return where(r > 1., 0., 93.3333*r**7 - 448.*r**6 + 840.*r**5 - 746.667*r**4 + 280.*r**3 - 18.6667*r)
def d_wendlandC6(r):     return where(r > 1., 0., 352.*r**10 - 2310.*r**9 + 6336.*r**8 - 9240.*r**7 + 7392.*r**6 - 2772.*r**5 + 264.*r**3 - 22.*r)
def d_Gaussian(r,h): return (-8.*sqrt(2.)/(PI**(3./2.) * h**5)) * r * exp(- 2.*r**2 / (h**2))

# Get the second derivative of the reduced kernel definitions for 3D kernels
def d2_cubic_spline(r):   return where(r > 1., 0., where(r < 0.5,
                                                         18.*r - 6.,
                                                         -6.*r + 6.) )

def d2_quartic_spline(r): return where(r > 1., 0., where(r < 0.2,
                                                         72.*r**2 - 4.8,
                                                         where(r < 0.6,
                                                               -48.*r**2 + 48.*r  - (48./5.),
                                                               12.*r**2 - 24.*r + 12.)))

def d2_quintic_spline(r): return where(r > 1., 0., where(r < 1./3.,
                                                         -200.*r**3 + 120.*r**2 - (40./9.),
                                                         where(r < 2./3.,
                                                               100.*r**3 - 180.*r**2 + 100.*r - (140./9.),
                                                               -20.*r**3 + 60.*r**2 - 60.*r + 20)))
def d2_wendlandC2(r):     return where(r > 1., 0., 80.*r**3 - 180.*r**2 + 120.*r - 20.)
def d2_wendlandC4(r):     return where(r > 1., 0., 653.3333*r**6 - 2688.*r**5 + 4200.*r**4 - 2986.667*r**3 + 840.*r**2 - 18.6667)
def d2_wendlandC6(r):     return where(r > 1., 0., 3520.*r**9 - 20790.*r**8 + 50688.*r**7 - 64680.*r**6 + 44352.*r**5 - 13860.*r**4 + 792.*r**2 - 22)
def d2_Gaussian(r,h): return (32*sqrt(2)/(PI**(3./2.)*h**7)) * r**2 * exp(-2.*r**2 / (h**2)) - 8.*sqrt(2.)/(PI**(3./2.) * h**5) * exp(- 2.*r**2 / (h**2))


# Derivative of kernel definitions (3D)
def dW_cubic_spline(r):   return C_cubic      * d_cubic_spline(r / H_cubic)     / H_cubic**4
def dW_quartic_spline(r): return C_quartic    * d_quartic_spline(r / H_quartic) / H_quartic**4
def dW_quintic_spline(r): return C_quintic    * d_quintic_spline(r / H_quintic) / H_quintic**4
def dW_WendlandC2(r):     return C_WendlandC2 * d_wendlandC2(r / H_WendlandC2)  / H_WendlandC2**4
def dW_WendlandC4(r):     return C_WendlandC4 * d_wendlandC4(r / H_WendlandC4)  / H_WendlandC4**4
def dW_WendlandC6(r):     return C_WendlandC6 * d_wendlandC6(r / H_WendlandC6)  / H_WendlandC6**4

# Second derivative of kernel definitions (3D)
def d2W_cubic_spline(r):   return C_cubic      * d2_cubic_spline(r / H_cubic)     / H_cubic**5
def d2W_quartic_spline(r): return C_quartic    * d2_quartic_spline(r / H_quartic) / H_quartic**5
def d2W_quintic_spline(r): return C_quintic    * d2_quintic_spline(r / H_quintic) / H_quintic**5
def d2W_WendlandC2(r):     return C_WendlandC2 * d2_wendlandC2(r / H_WendlandC2)  / H_WendlandC2**5
def d2W_WendlandC4(r):     return C_WendlandC4 * d2_wendlandC4(r / H_WendlandC4)  / H_WendlandC4**5
def d2W_WendlandC6(r):     return C_WendlandC6 * d2_wendlandC6(r / H_WendlandC6)  / H_WendlandC6**5


figure()
subplot(211)

plot([0, 2.5*h], [0., 0.], 'k--', linewidth=0.7)
plot(xx, d_Gaussian(xx, h), 'k-', linewidth=0.7, label="${\\rm Gaussian}$")
plot(xx, dW_cubic_spline(xx), 'b-', label="${\\rm Cubic~spline}$")
plot(xx, dW_quartic_spline(xx), 'c-', label="${\\rm Quartic~spline}$")
plot(xx, dW_quintic_spline(xx), 'g-', label="${\\rm Quintic~spline}$")
plot(xx, dW_WendlandC2(xx), 'r-', label="${\\rm Wendland~C2}$")
plot(xx, dW_WendlandC4(xx), 'm-', label="${\\rm Wendland~C4}$")
plot(xx, dW_WendlandC6(xx), 'y-', label="${\\rm Wendland~C6}$")

maxY = d_Gaussian(h/2, h)

# Show h
plot([h, h], [2*maxY, 0.1], 'k:', linewidth=0.5)

# Show sigma
plot([h/2, h/2], [2*maxY, 0.1], 'k:', linewidth=0.5)

# Show <x>
plot([dx, dx], [2*maxY, 0.1], 'k:', linewidth=0.5)


xlim(0., 2.5*h)
gca().xaxis.set_ticklabels([])
ylim(1.1*maxY, -0.1*maxY)
xlabel("$r$", labelpad=0)
ylabel("$\\partial W(r,h)/\\partial r$", labelpad=0.5)
legend(loc="lower right")



subplot(212)

maxY = d2_Gaussian(h,h)
plot([h/2, h/2], [-4*maxY, 1.4*maxY], 'k:', linewidth=0.5)
plot([h, h], [-4*maxY, 1.4*maxY], 'k:', linewidth=0.5)
plot([dx, dx], [-4*maxY, 1.4*maxY], 'k:', linewidth=0.5)
text(h/2, -3.*maxY, "$\\sigma\\equiv h/2$", rotation=90, backgroundcolor='w', ha='center', va='bottom')
text(h, -3.*maxY, "$h\\equiv\\eta\\langle x\\rangle = %.4f$"%h, rotation=90, backgroundcolor='w', ha='center', va='bottom')
text(dx, -3.*maxY, "$\\langle x\\rangle = %.1f$"%dx, rotation=90, backgroundcolor='w', ha='center', va='bottom')


plot([0, 2.5*h], [0., 0.], 'k--', linewidth=0.7)
plot(xx, d2_Gaussian(xx, h), 'k-', linewidth=0.7, label="${\\rm Gaussian}$")
plot(xx, d2W_cubic_spline(xx), 'b-', label="${\\rm Cubic~spline}$")
plot(xx, d2W_quartic_spline(xx), 'c-', label="${\\rm Quartic~spline}$")
plot(xx, d2W_quintic_spline(xx), 'g-', label="${\\rm Quintic~spline}$")
plot(xx, d2W_WendlandC2(xx), 'r-', label="${\\rm Wendland~C2}$")
plot(xx, d2W_WendlandC4(xx), 'm-', label="${\\rm Wendland~C4}$")
plot(xx, d2W_WendlandC6(xx), 'y-', label="${\\rm Wendland~C6}$")

xlim(0., 2.5*h)
ylim(-3.2*maxY, 1.3*maxY)
xlabel("$r$", labelpad=0)
ylabel("$\\partial^2 W(r,h)/\\partial r^2$", labelpad=0.5)


savefig("kernel_derivatives.pdf")
