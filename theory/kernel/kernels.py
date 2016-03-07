#!/usr/bin/env python2                                  
# -*- coding: utf-8 -*-            
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
    'axes.labelsize': 8,
    'axes.titlesize': 8,
    'font.size': 8,
    'legend.fontsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'xtick.major.pad': 2.5,
    'ytick.major.pad': 2.5,
    'text.usetex': True,
'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.12,
'figure.subplot.right'   : 0.99  ,
'figure.subplot.bottom'  : 0.08  ,
'figure.subplot.top'     : 0.99  ,
'figure.subplot.wspace'  : 0.  ,
'figure.subplot.hspace'  : 0.  ,
'lines.markersize' : 6,
'lines.linewidth' : 2.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


#Parameters
eta = 1.2349
h = 2.1

#Constants
PI = math.pi

#Cubic Spline
cubic_kernel_degree = 3
cubic_kernel_ivals = 2
cubic_kernel_gamma = 2.
cubic_kernel_ngb = 4.0 / 3.0 * PI * eta**3 * 6.0858
cubic_kernel_coeffs = array([[3./(4.*PI) , -3./(2.*PI), 0.,     1./PI],
                             [-1./(4.*PI),  3./(2.*PI), -3./PI, 2./PI],
                             [0.,           0.,         0.,     0.]])
def cubic_W(x):
    if size(x) == 1:
        x = array([x])
    ind = (minimum(x, cubic_kernel_ivals)).astype(int)
    coeffs = cubic_kernel_coeffs[ind,:]
    w = coeffs[:,0] * x + coeffs[:,1]
    for k in range(2, cubic_kernel_degree+1):
        w = x * w + coeffs[:,k]
    return w


#Quartic Spline
quartic_kernel_degree = 4
quartic_kernel_ivals = 3
quartic_kernel_gamma = 2.5
quartic_kernel_ngb = 4.0 / 3.0 * PI * eta**3 * 8.2293
quartic_kernel_coeffs = array([[3./(10.*PI) , 0.,           -3./(4.*PI) , 0.          ,  23./(32.*PI)],
                               [-1./(5.*PI) , 1./PI       , -3./(2.*PI) ,1./(4.*PI)   ,  11./(16.*PI)],
                               [1./(20.*PI) , -1./(2.*PI) , 15./(8.*PI) , -25./(8.*PI), 125./(64.*PI)],
                               [ 0. , 0.,           0.,         0.,     0.]])
def quartic_W(x):
    if size(x) == 1:
        x = array([x])
    ind = (minimum(x+0.5, quartic_kernel_ivals)).astype(int)
    coeffs = quartic_kernel_coeffs[ind,:]
    w = coeffs[:,0] * x + coeffs[:,1]
    for k in range(2, quartic_kernel_degree+1):
        w = x * w + coeffs[:,k]
    return w


# Wendland kernel
wendland2_kernel_degree = 5
wendland2_kernel_ivals = 1
wendland2_kernel_gamma = 2
wendland2_kernel_ngb = 4.0 / 3.0 * PI * eta**3 * 7.261825
wendland2_kernel_coeffs = 3.342253804929802 * array([[1./8, -30./32, 80./32, -80./32., 0., 1.],
                                                     [ 0. , 0.,  0.,   0., 0., 0.]]) / 8.

print wendland2_kernel_coeffs
def wendland2_W(x):
    if size(x) == 1:
        x = array([x])
    ind = (minimum(0.5*x, wendland2_kernel_ivals)).astype(int)
    coeffs = wendland2_kernel_coeffs[ind,:]
    w = coeffs[:,0] * x + coeffs[:,1]
    for k in range(2, wendland2_kernel_degree+1):
        w = x * w + coeffs[:,k]
    return w 

#def wendland2_W(x):
#    if size(x) == 1:
#        x = array([x])
#    x /= 1.936492
#    x[x>1.] = 1.
#    oneminusu = 1.-x
#    oneminusu4 = oneminusu * oneminusu * oneminusu * oneminusu
#    return 3.342253804929802 * oneminusu4 * (1. + 4.*x) / h**3


#Find H
r = arange(0, 3.5*h, 1./1000.)
xi = r/h
cubic_Ws = cubic_W(xi)
quartic_Ws = quartic_W(xi)
wendland2_Ws = wendland2_W(xi)
for j in range(size(r)):
    if cubic_Ws[j] == 0:
        cubic_H = r[j]
        break
for j in range(size(r)):
    if quartic_Ws[j] == 0:
        quartic_H = r[j]
        break
for j in range(size(r)):
    if wendland2_Ws[j] == 0:
        wendland2_H = r[j]
        break

    
print "H=", cubic_H
print "H=", quartic_H
print "H=", wendland2_H


# Compute sigma -----------------------------------------
cubic_norm = 4.*PI*integrate.quad(lambda x: x**2*cubic_W(x), 0, cubic_H)[0]
quartic_norm = 4.*PI*integrate.quad(lambda x: x**2*quartic_W(x), 0, quartic_H)[0]
wendland2_norm = 4.*PI*integrate.quad(lambda x: x**2*wendland2_W(x), 0, wendland2_H)[0]

print cubic_norm
print quartic_norm
print wendland2_norm


# Plot kernels ------------------------------------------
r = arange(0, 3.5*h, 1./100.)
xi = r/h

cubic_Ws = cubic_W(xi)
quartic_Ws = quartic_W(xi)
wendland2_Ws = wendland2_W(xi)



figure()

text(h-0.1, cubic_Ws[0]/20., "h", ha="right",va="center")
arrow(h, cubic_Ws[0]/10., 0., -cubic_Ws[0]/10., fc='k', ec='k', length_includes_head=True, head_length=cubic_Ws[0]/30., head_width=0.1)


plot(r,cubic_Ws, 'b-' ,label="Cubic")
plot(r, quartic_Ws, 'r-', label="Quartic")
plot(r, wendland2_Ws, 'g-', label="Wendland C2")

text(cubic_H-0.1, cubic_Ws[0]/20., "H", ha="right",va="center", color='b')
arrow(cubic_H, cubic_Ws[0]/10., 0., -cubic_Ws[0]/10., fc='b', ec='b', length_includes_head=True, head_length=cubic_Ws[0]/30., head_width=0.1)

text(quartic_H-0.1, cubic_Ws[0]/20., "H", ha="right",va="center", color='r')
arrow(quartic_H, cubic_Ws[0]/10., 0., -cubic_Ws[0]/10., fc='r', ec='r', length_includes_head=True, head_length=quartic_Ws[0]/30., head_width=0.1)

text(wendland2_H-0.1, cubic_Ws[0]/20., "H", ha="right",va="center", color='r')
arrow(wendland2_H, cubic_Ws[0]/10., 0., -cubic_Ws[0]/10., fc='g', ec='g', length_includes_head=True, head_length=wendland2_Ws[0]/30., head_width=0.1)


xlabel("r", labelpad=0)
ylabel("W(r,h)", labelpad=0)
legend(loc="upper right")
savefig("kernel.pdf")

