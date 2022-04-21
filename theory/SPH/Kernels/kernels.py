###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from matplotlib.font_manager import FontProperties
import numpy


params = {
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "text.usetex": True,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.17,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.08,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})


# params = {
#     'axes.labelsize': 10,
#     'axes.titlesize': 8,
#     'font.size': 10,
#     'legend.fontsize': 9,
#     'xtick.labelsize': 10,
#     'ytick.labelsize': 10,
#     'xtick.major.pad': 2.5,
#     'ytick.major.pad': 2.5,
#     'text.usetex': True,
# 'figure.figsize' : (4.15,4.15),
# 'figure.subplot.left'    : 0.14,
# 'figure.subplot.right'   : 0.99,
# 'figure.subplot.bottom'  : 0.06,
# 'figure.subplot.top'     : 0.99,
# 'figure.subplot.wspace'  : 0.  ,
# 'figure.subplot.hspace'  : 0.  ,
# 'lines.markersize' : 6,
# 'lines.linewidth' : 1.5,
# 'text.latex.unicode': True
# }
# rcParams.update(params)
# rc('font',**{'family':'sans-serif','sans-serif':['Times']})


# Parameters
eta = 1.2348422195325  # Resolution (Gives 48 neighbours for a cubic spline kernel)
dx = 1.5  # 4 #2.7162  # Mean inter-particle separation

# Constants
PI = math.pi

# Compute expected moothing length
h = eta * dx

# Get kernel support (Dehnen & Aly 2012, table 1) for 3D kernels
H_cubic = 1.825742 * h
H_quartic = 2.018932 * h
H_quintic = 2.195775 * h
H_WendlandC2 = 1.936492 * h
H_WendlandC4 = 2.207940 * h
H_WendlandC6 = 2.449490 * h

# Get the number of neighbours within kernel support:
N_H_cubic = 4.0 / 3.0 * PI * H_cubic ** 3 / (dx) ** 3
N_H_quartic = 4.0 / 3.0 * PI * H_quartic ** 3 / (dx) ** 3
N_H_quintic = 4.0 / 3.0 * PI * H_quintic ** 3 / (dx) ** 3
N_H_WendlandC2 = 4.0 / 3.0 * PI * H_WendlandC2 ** 3 / (dx) ** 3
N_H_WendlandC4 = 4.0 / 3.0 * PI * H_WendlandC4 ** 3 / (dx) ** 3
N_H_WendlandC6 = 4.0 / 3.0 * PI * H_WendlandC6 ** 3 / (dx) ** 3


print(
    "Smoothing length: h =",
    h,
    "Cubic spline kernel support size:   H =",
    H_cubic,
    "Number of neighbours N_H =",
    N_H_cubic,
)
print(
    "Smoothing length: h =",
    h,
    "Quartic spline kernel support size: H =",
    H_quartic,
    "Number of neighbours N_H =",
    N_H_quartic,
)
print(
    "Smoothing length: h =",
    h,
    "Quintic spline kernel support size: H =",
    H_quintic,
    "Number of neighbours N_H =",
    N_H_quintic,
)
print(
    "Smoothing length: h =",
    h,
    "Wendland C2 kernel support size:    H =",
    H_WendlandC2,
    "Number of neighbours N_H =",
    N_H_WendlandC2,
)
print(
    "Smoothing length: h =",
    h,
    "Wendland C4 kernel support size:    H =",
    H_WendlandC4,
    "Number of neighbours N_H =",
    N_H_WendlandC4,
)
print(
    "Smoothing length: h =",
    h,
    "Wendland C6 kernel support size:    H =",
    H_WendlandC6,
    "Number of neighbours N_H =",
    N_H_WendlandC6,
)

# Get kernel constants (Dehen & Aly 2012, table 1) for 3D kernel
C_cubic = 16.0 / PI
C_quartic = 5 ** 6 / (512 * PI)
C_quintic = 3 ** 7 / (40 * PI)
C_WendlandC2 = 21.0 / (2 * PI)
C_WendlandC4 = 495.0 / (32 * PI)
C_WendlandC6 = 1365.0 / (64 * PI)

# Get the reduced kernel definitions (Dehen & Aly 2012, table 1) for 3D kernel
# def plus(u) : return maximum(0., u)
def cubic_spline(r):
    return where(
        r > 1.0,
        0.0,
        where(
            r < 0.5,
            3.0 * r ** 3 - 3.0 * r ** 2 + 0.5,
            -r ** 3 + 3.0 * r ** 2 - 3.0 * r + 1.0,
        ),
    )


# return plus(1. - r)**3 - 4.*plus(1./2. - r)**3
def quartic_spline(r):
    return where(
        r > 1.0,
        0.0,
        where(
            r < 0.2,
            6.0 * r ** 4 - 2.4 * r ** 2 + 46.0 / 125.0,
            where(
                r < 0.6,
                -4.0 * r ** 4
                + 8.0 * r ** 3
                - (24.0 / 5.0) * r ** 2
                + (8.0 / 25.0) * r
                + 44.0 / 125.0,
                1.0 * r ** 4 - 4.0 * r ** 3 + 6.0 * r ** 2 - 4.0 * r + 1.0,
            ),
        ),
    )


# return plus(1. - r)**4 - 5.*plus(3./5. - r)**4 + 10.*plus(1./5. - r)**4
def quintic_spline(r):
    return where(
        r > 1.0,
        0.0,
        where(
            r < 1.0 / 3.0,
            -10.0 * r ** 5 + 10.0 * r ** 4 - (20.0 / 9.0) * r ** 2 + (22.0 / 81.0),
            where(
                r < 2.0 / 3.0,
                5.0 * r ** 5
                - 15.0 * r ** 4
                + (50.0 / 3.0) * r ** 3
                - (70.0 / 9.0) * r ** 2
                + (25.0 / 27.0) * r
                + (17.0 / 81.0),
                -1.0 * r ** 5
                + 5.0 * r ** 4
                - 10.0 * r ** 3
                + 10.0 * r ** 2
                - 5.0 * r
                + 1.0,
            ),
        ),
    )


# return plus(1. - r)**5 - 6.*plus(2./3. - r)**5 + 15.*plus(1./3. - r)**5
def wendlandC2(r):
    return where(
        r > 1.0, 0.0, 4.0 * r ** 5 - 15.0 * r ** 4 + 20.0 * r ** 3 - 10 * r ** 2 + 1.0
    )


def wendlandC4(r):
    return where(
        r > 1.0,
        0.0,
        (35.0 / 3.0) * r ** 8
        - 64.0 * r ** 7
        + 140.0 * r ** 6
        - (448.0 / 3.0) * r ** 5
        + 70.0 * r ** 4
        - (28.0 / 3.0) * r ** 2
        + 1.0,
    )


def wendlandC6(r):
    return where(
        r > 1.0,
        0.0,
        32.0 * r ** 11
        - 231.0 * r ** 10
        + 704.0 * r ** 9
        - 1155.0 * r ** 8
        + 1056.0 * r ** 7
        - 462.0 * r ** 6
        + 66.0 * r ** 4
        - 11.0 * r ** 2
        + 1.0,
    )


def Gaussian(r, h):
    return (1.0 / (0.5 * pi * h ** 2) ** (3.0 / 2.0)) * exp(-2.0 * r ** 2 / (h ** 2))


# Kernel definitions (3D)
def W_cubic_spline(r):
    return C_cubic * cubic_spline(r / H_cubic) / H_cubic ** 3


def W_quartic_spline(r):
    return C_quartic * quartic_spline(r / H_quartic) / H_quartic ** 3


def W_quintic_spline(r):
    return C_quintic * quintic_spline(r / H_quintic) / H_quintic ** 3


def W_WendlandC2(r):
    return C_WendlandC2 * wendlandC2(r / H_WendlandC2) / H_WendlandC2 ** 3


def W_WendlandC4(r):
    return C_WendlandC4 * wendlandC4(r / H_WendlandC4) / H_WendlandC4 ** 3


def W_WendlandC6(r):
    return C_WendlandC6 * wendlandC6(r / H_WendlandC6) / H_WendlandC6 ** 3


# PLOT STUFF
figure()
subplot(211)
xx = linspace(0.0, 5 * h, 100)
maxY = 1.2 * Gaussian(0, h)

# Plot the kernels
plot(
    xx,
    Gaussian(xx, h),
    "k-",
    linewidth=0.7,
    label="${\\rm %14s~ H=\\infty}$" % ("Gaussian~~~~~~"),
)
plot(
    xx,
    W_cubic_spline(xx),
    "b-",
    label="${\\rm %14s~ H=%4.3f}$" % ("Cubic~spline~~", H_cubic),
    lw=1.5,
)
plot(
    xx,
    W_quartic_spline(xx),
    "c-",
    label="${\\rm %14s~ H=%4.3f}$" % ("Quartic~spline", H_quartic),
    lw=1.5,
)
plot(
    xx,
    W_quintic_spline(xx),
    "g-",
    label="${\\rm %14s~ H=%4.3f}$" % ("Quintic~spline", H_quintic),
    lw=1.5,
)
plot(
    xx,
    W_WendlandC2(xx),
    "r-",
    label="${\\rm %14s~ H=%4.3f}$" % ("Wendland~C2~", H_WendlandC2),
    lw=1.5,
)
plot(
    xx,
    W_WendlandC4(xx),
    "m-",
    label="${\\rm %14s~ H=%4.3f}$" % ("Wendland~C4~", H_WendlandC4),
    lw=1.5,
)
plot(
    xx,
    W_WendlandC6(xx),
    "y-",
    label="${\\rm %14s~ H=%4.3f}$" % ("Wendland~C6~", H_WendlandC6),
    lw=1.5,
)

# Indicate the position of H
arrow(
    H_cubic,
    0.12 * maxY,
    0.0,
    -0.12 * maxY * 0.9,
    fc="b",
    ec="b",
    length_includes_head=True,
    head_width=0.03,
    head_length=0.12 * maxY * 0.3,
)
arrow(
    H_quartic,
    0.12 * maxY,
    0.0,
    -0.12 * maxY * 0.9,
    fc="c",
    ec="c",
    length_includes_head=True,
    head_width=0.03,
    head_length=0.12 * maxY * 0.3,
)
arrow(
    H_quintic,
    0.12 * maxY,
    0.0,
    -0.12 * maxY * 0.9,
    fc="g",
    ec="g",
    length_includes_head=True,
    head_width=0.03,
    head_length=0.12 * maxY * 0.3,
)
arrow(
    H_WendlandC2,
    0.12 * maxY,
    0.0,
    -0.12 * maxY * 0.9,
    fc="r",
    ec="r",
    length_includes_head=True,
    head_width=0.03,
    head_length=0.12 * maxY * 0.3,
)
arrow(
    H_WendlandC4,
    0.12 * maxY,
    0.0,
    -0.12 * maxY * 0.9,
    fc="m",
    ec="m",
    length_includes_head=True,
    head_width=0.03,
    head_length=0.12 * maxY * 0.3,
)
arrow(
    H_WendlandC6,
    0.12 * maxY,
    0.0,
    -0.12 * maxY * 0.9,
    fc="y",
    ec="y",
    length_includes_head=True,
    head_width=0.03,
    head_length=0.12 * maxY * 0.3,
)

# Show h
plot([h, h], [0.0, 0.63 * maxY], "k:", linewidth=0.5)
text(
    h,
    maxY * 0.65,
    "$h\\equiv\\eta\\langle x\\rangle$",
    rotation=90,
    ha="center",
    va="bottom",
    fontsize=9,
)

# Show sigma
plot([h / 2, h / 2], [0.0, 0.63 * maxY], "k:", linewidth=0.5)
text(
    h / 2,
    maxY * 0.65,
    "$\\sigma\\equiv h/2$",
    rotation=90,
    ha="center",
    va="bottom",
    fontsize=9,
)

# Show <x>
plot([dx, dx], [0.0, 0.63 * maxY], "k:", linewidth=0.5)
text(
    dx,
    maxY * 0.65,
    "$\\langle x\\rangle = %.1f$" % dx,
    rotation=90,
    ha="center",
    va="bottom",
    fontsize=9,
)

xlim(0.0, 2.5 * h)
ylim(0.0, maxY)
gca().xaxis.set_ticklabels([])
ylabel("$W(r,h)$", labelpad=1.5)
legend(
    loc="upper right", frameon=False, handletextpad=0.1, handlelength=1.2, fontsize=8
)


# Same but now in log space
subplot(212, yscale="log")
plot(xx, Gaussian(xx, h), "k-", linewidth=0.7, label="${\\rm Gaussian}$", lw=1.5)
plot(xx, W_cubic_spline(xx), "b-", label="${\\rm Cubic~spline}$", lw=1.5)
plot(xx, W_quartic_spline(xx), "c-", label="${\\rm Quartic~spline}$", lw=1.5)
plot(xx, W_quintic_spline(xx), "g-", label="${\\rm Quintic~spline}$", lw=1.5)
plot(xx, W_WendlandC2(xx), "r-", label="${\\rm Wendland~C2}$", lw=1.5)
plot(xx, W_WendlandC4(xx), "m-", label="${\\rm Wendland~C4}$", lw=1.5)
plot(xx, W_WendlandC6(xx), "y-", label="${\\rm Wendland~C6}$", lw=1.5)

# Show h
plot([h, h], [0.0, 1.0], "k:", linewidth=0.5)

# Show sigma
plot([h / 2, h / 2], [0.0, 1.0], "k:", linewidth=0.5)

# Show <x>
plot([dx, dx], [0.0, 1.0], "k:", linewidth=0.5)


# Show plot properties
text(
    h / 5.0,
    2e-3,
    "$\\langle x \\rangle = %3.1f$" % (dx),
    va="top",
    backgroundcolor="w",
    fontsize=9,
)
text(
    h / 5.0 + 0.06,
    4e-4,
    "$\\eta = %5.4f$" % (eta),
    va="top",
    backgroundcolor="w",
    fontsize=9,
)
text(
    h / 5.0 + 0.06,
    8e-5,
    "$h = \\eta\\langle x \\rangle = %5.4f$" % (h),
    va="top",
    backgroundcolor="w",
    fontsize=9,
)

# Show number of neighbours
text(1.75 * h, 2e-1 / 2.9 ** 0, "$N_{\\rm ngb}=\\infty$", fontsize=10)
text(
    1.75 * h,
    2e-1 / 2.9 ** 1,
    "$N_{\\rm ngb}=%3.1f$" % (N_H_cubic),
    color="b",
    fontsize=9,
)
text(
    1.75 * h,
    2e-1 / 2.9 ** 2,
    "$N_{\\rm ngb}=%3.1f$" % (N_H_quartic),
    color="c",
    fontsize=9,
)
text(
    1.75 * h,
    2e-1 / 2.9 ** 3,
    "$N_{\\rm ngb}=%3.1f$" % (N_H_quintic),
    color="g",
    fontsize=9,
)
text(
    1.75 * h,
    2e-1 / 2.9 ** 4,
    "$N_{\\rm ngb}=%3.1f$" % (N_H_WendlandC2),
    color="r",
    fontsize=9,
)
text(
    1.75 * h,
    2e-1 / 2.9 ** 5,
    "$N_{\\rm ngb}=%3.1f$" % (N_H_WendlandC4),
    color="m",
    fontsize=9,
)
text(
    1.75 * h,
    2e-1 / 2.9 ** 6,
    "$N_{\\rm ngb}=%3.0f$" % (N_H_WendlandC6),
    color="y",
    fontsize=9,
)

xlim(0.0, 2.5 * h)
ylim(1e-5, 0.7)
xlabel("$r$", labelpad=-1)
ylabel("$W(r,h)$", labelpad=0.5)

savefig("kernels.pdf")


################################
# Now, let's work on derivatives
################################

# Get the derivative of the reduced kernel definitions for 3D kernels
def d_cubic_spline(r):
    return where(
        r > 1.0,
        0.0,
        where(r < 0.5, 9.0 * r ** 2 - 6.0 * r, -3.0 * r ** 2 + 6.0 * r - 3.0),
    )


def d_quartic_spline(r):
    return where(
        r > 1.0,
        0.0,
        where(
            r < 0.2,
            24.0 * r ** 3 - 4.8 * r,
            where(
                r < 0.6,
                -16.0 * r ** 3 + 24.0 * r ** 2 - (48.0 / 5.0) * r + (8.0 / 25.0),
                4.0 * r ** 3 - 12.0 * r ** 2 + 12.0 * r - 4.0,
            ),
        ),
    )


def d_quintic_spline(r):
    return where(
        r > 1.0,
        0.0,
        where(
            r < 1.0 / 3.0,
            -50.0 * r ** 4 + 40.0 * r ** 3 - (40.0 / 9.0) * r,
            where(
                r < 2.0 / 3.0,
                25.0 * r ** 4
                - 60.0 * r ** 3
                + 50.0 * r ** 2
                - (140.0 / 9.0) * r
                + (25.0 / 27.0),
                -5.0 * r ** 4 + 20.0 * r ** 3 - 30.0 * r ** 2 + 20.0 * r - 5.0,
            ),
        ),
    )


def d_wendlandC2(r):
    return where(r > 1.0, 0.0, 20.0 * r ** 4 - 60.0 * r ** 3 + 60.0 * r ** 2 - 20.0 * r)


def d_wendlandC4(r):
    return where(
        r > 1.0,
        0.0,
        93.3333 * r ** 7
        - 448.0 * r ** 6
        + 840.0 * r ** 5
        - 746.667 * r ** 4
        + 280.0 * r ** 3
        - 18.6667 * r,
    )


def d_wendlandC6(r):
    return where(
        r > 1.0,
        0.0,
        352.0 * r ** 10
        - 2310.0 * r ** 9
        + 6336.0 * r ** 8
        - 9240.0 * r ** 7
        + 7392.0 * r ** 6
        - 2772.0 * r ** 5
        + 264.0 * r ** 3
        - 22.0 * r,
    )


def d_Gaussian(r, h):
    return (
        (-8.0 * sqrt(2.0) / (PI ** (3.0 / 2.0) * h ** 5))
        * r
        * exp(-2.0 * r ** 2 / (h ** 2))
    )


def dh_Gaussian(r, h):
    return -(3 * Gaussian(r, h) + (r / h) * d_Gaussian(r, h))


# Get the second derivative of the reduced kernel definitions for 3D kernels
# def d2_cubic_spline(r):   return where(r > 1., 0., where(r < 0.5,
#                                                          18.*r - 6.,
#                                                          -6.*r + 6.) )

# def d2_quartic_spline(r): return where(r > 1., 0., where(r < 0.2,
#                                                          72.*r**2 - 4.8,
#                                                          where(r < 0.6,
#                                                                -48.*r**2 + 48.*r  - (48./5.),
#                                                                12.*r**2 - 24.*r + 12.)))

# def d2_quintic_spline(r): return where(r > 1., 0., where(r < 1./3.,
#                                                          -200.*r**3 + 120.*r**2 - (40./9.),
#                                                          where(r < 2./3.,
#                                                                100.*r**3 - 180.*r**2 + 100.*r - (140./9.),
#                                                                -20.*r**3 + 60.*r**2 - 60.*r + 20)))
# def d2_wendlandC2(r):     return where(r > 1., 0., 80.*r**3 - 180.*r**2 + 120.*r - 20.)
# def d2_wendlandC4(r):     return where(r > 1., 0., 653.3333*r**6 - 2688.*r**5 + 4200.*r**4 - 2986.667*r**3 + 840.*r**2 - 18.6667)
# def d2_wendlandC6(r):     return where(r > 1., 0., 3520.*r**9 - 20790.*r**8 + 50688.*r**7 - 64680.*r**6 + 44352.*r**5 - 13860.*r**4 + 792.*r**2 - 22)
# def d2_Gaussian(r,h): return (32*sqrt(2)/(PI**(3./2.)*h**7)) * r**2 * exp(-2.*r**2 / (h**2)) - 8.*sqrt(2.)/(PI**(3./2.) * h**5) * exp(- 2.*r**2 / (h**2))


# Derivative of kernel definitions (3D)
def dWdx_cubic_spline(r):
    return C_cubic * d_cubic_spline(r / H_cubic) / H_cubic ** 4


def dWdx_quartic_spline(r):
    return C_quartic * d_quartic_spline(r / H_quartic) / H_quartic ** 4


def dWdx_quintic_spline(r):
    return C_quintic * d_quintic_spline(r / H_quintic) / H_quintic ** 4


def dWdx_WendlandC2(r):
    return C_WendlandC2 * d_wendlandC2(r / H_WendlandC2) / H_WendlandC2 ** 4


def dWdx_WendlandC4(r):
    return C_WendlandC4 * d_wendlandC4(r / H_WendlandC4) / H_WendlandC4 ** 4


def dWdx_WendlandC6(r):
    return C_WendlandC6 * d_wendlandC6(r / H_WendlandC6) / H_WendlandC6 ** 4


# Derivative of kernel definitions (3D)
def dWdh_cubic_spline(r):
    return (
        3.0 * cubic_spline(r / H_cubic) + (r / H_cubic) * d_cubic_spline(r / H_cubic)
    ) * (-C_cubic / H_cubic ** 4)


def dWdh_quartic_spline(r):
    return (
        3.0 * quartic_spline(r / H_quartic)
        + (r / H_quartic) * d_quartic_spline(r / H_quartic)
    ) * (-C_quartic / H_quartic ** 4)


def dWdh_quintic_spline(r):
    return (
        3.0 * quintic_spline(r / H_quintic)
        + (r / H_quintic) * d_quintic_spline(r / H_quintic)
    ) * (-C_quintic / H_quintic ** 4)


def dWdh_WendlandC2(r):
    return (
        3.0 * wendlandC2(r / H_WendlandC2)
        + (r / H_WendlandC2) * d_wendlandC2(r / H_WendlandC2)
    ) * (-C_WendlandC2 / H_WendlandC2 ** 4)


def dWdh_WendlandC4(r):
    return (
        3.0 * wendlandC4(r / H_WendlandC4)
        + (r / H_WendlandC4) * d_wendlandC4(r / H_WendlandC4)
    ) * (-C_WendlandC4 / H_WendlandC4 ** 4)


def dWdh_WendlandC6(r):
    return (
        3.0 * wendlandC6(r / H_WendlandC6)
        + (r / H_WendlandC6) * d_wendlandC6(r / H_WendlandC6)
    ) * (-C_WendlandC6 / H_WendlandC6 ** 4)


# Second derivative of kernel definitions (3D)
# def d2W_cubic_spline(r):   return C_cubic      * d2_cubic_spline(r / H_cubic)     / H_cubic**5
# def d2W_quartic_spline(r): return C_quartic    * d2_quartic_spline(r / H_quartic) / H_quartic**5
# def d2W_quintic_spline(r): return C_quintic    * d2_quintic_spline(r / H_quintic) / H_quintic**5
# def d2W_WendlandC2(r):     return C_WendlandC2 * d2_wendlandC2(r / H_WendlandC2)  / H_WendlandC2**5
# def d2W_WendlandC4(r):     return C_WendlandC4 * d2_wendlandC4(r / H_WendlandC4)  / H_WendlandC4**5
# def d2W_WendlandC6(r):     return C_WendlandC6 * d2_wendlandC6(r / H_WendlandC6)  / H_WendlandC6**5


figure()
subplot(211)

plot([0, 2.5 * h], [0.0, 0.0], "k--", linewidth=0.7)
plot(xx, d_Gaussian(xx, h), "k-", linewidth=0.7, label="${\\rm Gaussian}$", lw=1.5)
plot(xx, dWdx_cubic_spline(xx), "b-", label="${\\rm Cubic~spline}$", lw=1.5)
plot(xx, dWdx_quartic_spline(xx), "c-", label="${\\rm Quartic~spline}$", lw=1.5)
plot(xx, dWdx_quintic_spline(xx), "g-", label="${\\rm Quintic~spline}$", lw=1.5)
plot(xx, dWdx_WendlandC2(xx), "r-", label="${\\rm Wendland~C2}$", lw=1.5)
plot(xx, dWdx_WendlandC4(xx), "m-", label="${\\rm Wendland~C4}$", lw=1.5)
plot(xx, dWdx_WendlandC6(xx), "y-", label="${\\rm Wendland~C6}$", lw=1.5)

maxY = d_Gaussian(h / 2, h)

# Show h
plot([h, h], [2 * maxY, 0.1], "k:", linewidth=0.5)

# Show sigma
plot([h / 2, h / 2], [2 * maxY, 0.1], "k:", linewidth=0.5)

# Show <x>
plot([dx, dx], [2 * maxY, 0.1], "k:", linewidth=0.5)


xlim(0.0, 2.5 * h)
gca().xaxis.set_ticklabels([])
ylim(1.3 * maxY, -0.1 * maxY)
xlabel("$r$", labelpad=0)
ylabel("$\\partial W(r,h)/\\partial r$", labelpad=0.5)
legend(
    loc="lower right", frameon=False, handletextpad=0.1, handlelength=1.2, fontsize=8
)


subplot(212)
plot([h / 2, h / 2], [0.77 * maxY, -0.15 * maxY], "k:", linewidth=0.5)
plot([h, h], [0.77 * maxY, -0.15 * maxY], "k:", linewidth=0.5)
plot([dx, dx], [0.77 * maxY, -0.15 * maxY], "k:", linewidth=0.5)
text(
    h / 2,
    1.25 * maxY,
    "$\\sigma\\equiv h/2$",
    rotation=90,
    ha="center",
    va="bottom",
    fontsize=9,
)
text(
    h,
    1.25 * maxY,
    "$h\\equiv\\eta\\langle x\\rangle$",
    rotation=90,
    ha="center",
    va="bottom",
    fontsize=9,
)
text(
    dx,
    1.25 * maxY,
    "$\\langle x\\rangle = %.1f$" % dx,
    rotation=90,
    ha="center",
    va="bottom",
    fontsize=9,
)


plot([0, 2.5 * h], [0.0, 0.0], "k--", linewidth=0.7)
# plot(xx, dh_Gaussian(xx, h), 'k-', linewidth=0.7, label="${\\rm Gaussian}$")
plot(xx, dWdh_cubic_spline(xx), "b-", label="${\\rm Cubic~spline}$", lw=1.5)
plot(xx, dWdh_quartic_spline(xx), "c-", label="${\\rm Quartic~spline}$", lw=1.5)
plot(xx, dWdh_quintic_spline(xx), "g-", label="${\\rm Quintic~spline}$", lw=1.5)
plot(xx, dWdh_WendlandC2(xx), "r-", label="${\\rm Wendland~C2}$", lw=1.5)
plot(xx, dWdh_WendlandC4(xx), "m-", label="${\\rm Wendland~C4}$", lw=1.5)
plot(xx, dWdh_WendlandC6(xx), "y-", label="${\\rm Wendland~C6}$", lw=1.5)

xlim(0.0, 2.5 * h)
ylim(1.3 * maxY, -0.15 * maxY)
xlabel("$r$", labelpad=-1)
ylabel("$\\partial W(r,h)/\\partial h$", labelpad=0.5)


savefig("kernel_derivatives.pdf")
