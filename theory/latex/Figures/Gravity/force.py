import matplotlib
matplotlib.use("Agg")
from pylab import *
# tex stuff
#rc('font',**{'family':'serif','serif':['Palatino']})
params = {'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14,
          'text.usetex': True,

          'figure.figsize' : (9,6),
          'figure.subplot.left'    : 0.08,  # the left side of the subplots of the figure
          'figure.subplot.right'   : 0.95  ,  # the right side of the subplots of the figure
          'figure.subplot.bottom'  : 0.1  ,  # the bottom of the subplots of the figure
          'figure.subplot.top'     : 0.93  ,  # the top of the subplots of the figure
          'figure.subplot.wspace'  : 0.2  ,  # the amount of width reserved for blank space between subplots
          'figure.subplot.hspace'  : 0.2  ,  # the amount of height reserved for white space between subplots

          'lines.markersize' : 4,

          #'axes.formatter.limits' : (, 0),

          'text.latex.unicode': True
        }
rcParams.update(params)
rc('font', family='serif')
import sys
import math
import os



epsilon = 1.
h = 2.8*epsilon
r_s = 20.
r_c = 4.5*r_s

r=arange(0.1, 400, 1/1000.)

numElem = size(r)
f = zeros(numElem)
fac = zeros(numElem)
kernel = zeros(numElem)

for i in range(numElem):
    if ( r[i] >= r_c ):
        f[i] = 0.
    elif ( r[i] >= h ):
        f[i] = (1. / r[i]**3)
    else:
        u = r[i] / h
        if( u < 0.5 ):
            f[i] = (10.666666666667 + u * u * (32.0 * u - 38.4)) / h**3 
        else:
            f[i] = (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u)) / h**3

    fac[i] = math.erfc( r[i] / ( 2. * r_s ) ) + ( r[i] / ( r_s * sqrt( math.pi ) ) ) * exp( - r[i] * r[i] / ( 4 * r_s * r_s ) ) 
    f[i] = f[i]*fac[i]

for i in range(numElem):
    u = r[i] / h
    if( u < 0.5 ):
        kernel[i] = (10.666666666667 + u * u * (32.0 * u - 38.4)) / h**3 
    else:
        kernel[i] = (21.333333333333 - 48.0 * u + 38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u)) / h**3

figure()
loglog(r/h, 1/r**3, 'b-', label="Newton's force")
loglog(r/h, kernel, 'r-', label="Softening kernel")
loglog(r/h, fac/r**3, 'g-', label="Tree force",)
loglog(r/h, f, 'k--', label="Total short-range force", linewidth=3)


plot([1, 1], [1e-10, 1e10], 'k--')
text(0.85, 7e-9, "$h=2.8\epsilon$", rotation="vertical", fontsize=14, va="bottom")

plot([epsilon/h, epsilon/h], [1e-10, 1e10], 'k--')
text(0.85*epsilon/h, 7e-9, "$\epsilon$", rotation="vertical", fontsize=14, va="bottom")

plot([r_s/h, r_s/h], [1e-10, 1e10], 'k--')
text(0.85*r_s/h, 7e-9, "$r_s$", rotation="vertical", fontsize=14, va="bottom")

plot([r_c/h, r_c/h], [1e-10, 1e10], 'k--')
text(0.85*r_c/h, 7e-9, "$r_c$", rotation="vertical", fontsize=14, va="bottom")


legend(loc="upper right")


grid()

xlim(1e-1, 200)
xlabel("r/h")

ylim(4e-10, 1e2)
ylabel("Normalized acceleration")


savefig("force.png")
