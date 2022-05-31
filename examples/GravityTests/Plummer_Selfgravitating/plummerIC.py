################################################################################
# Copyright (c) 2022 Patrick Hirling (patrick.hirling@epfl.ch)
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
################################################################################

import numpy as np
import scipy.special as sci
from scipy.optimize import minimize
import time
import argparse
from swiftsimio import Writer
import unyt

parser = argparse.ArgumentParser(
    description="Generate discrete realization of Anisotropic Plummer Model"
)

# Parse Main Arguments
parser.add_argument("-a", type=float, default=0.05, help="Softening Length (in kpc)")
parser.add_argument(
    "-M", type=float, default=1.0e-5, help="Total Mass of Model (in 10^10 solar masses)"
)
parser.add_argument(
    "-q", type=float, default=1.0, help="Anisotropy parameter (between -inf and +2)"
)
parser.add_argument(
    "-N", type=int, default=1e4, help="Number of particles in the simulation"
)
parser.add_argument(
    "-bound", type=float, default=2.0, help="Upper bound of radius sampling"
)
parser.add_argument("-boxsize", type=float, default=4.0, help="Box Size of simulation")
parser.add_argument(
    "--noparallel",
    action="store_false",
    help="Do not use multiprocessing library (default: false)",
)
parser.add_argument(
    "-nbthreads",
    type=int,
    default=6,
    help="Number of threads to use for multiprocessing (default: 6)",
)

args = parser.parse_args()
if args.q > 2.0:
    raise ValueError(
        "Anisotropy parameter must be smaller than 2.0 (q = {:.2f} given)".format(
            args.q
        )
    )

#### Parameters
# Plummer Model
G = 4.299581e04  # Gravitational constant [kpc / 10^10 M_s * (kms)^2]
a = args.a  # Plummer softening length [kpc]
M = args.M  # Total Mass [10^10 M_s]
q = args.q  # Anisotropy Parameter (-inf,2]

# IC File
N = int(args.N)  # Number of Particles
fname = "plummer.hdf5"  # Name of the ic file (dont forget .hdf5)
pickle_ics = 0  # Optional: Pickle ics (pos,vel,mass) for future use
periodic_bdry = 0

# Parallelism for velocity sampling
use_parallel = args.noparallel
nb_threads = int(args.nbthreads)  # Set to None to use os.cpucount()

### Print parameters
print(
    "\n############################################################################################# \n"
)
print("Generating Plummer ICs with the following parameters::\n")
print(
    "Softening length:                                  a = " + "{:.3e} kpc".format(a)
)
print(
    "Total Mass:                                        M = "
    + "{:.3e} Solar Masses".format(M * 1e10)
)
print("Anisotropy parameter:                              q = " + str(q))
print("Output file:                                       " + fname)
print("Number of particles:                               N = " + "{:.1e}".format(N))
print(
    "\n############################################################################################# \n"
)

# Estimate good softening length (density at origin, sphere s.t. contains on avg 1 particle)
epsilon0 = a * N ** (-1.0 / 3.0)
print(
    "Recommended Softening length (times conventional factor): "
    + "{:.4e}".format(epsilon0)
    + " kpc"
)

### Generate Positions
# Inverse transform sampling from analytical inverse of M(<r)
def r_of_m(m, a, M):
    return a * ((M / m) ** (2.0 / 3.0) - 1) ** (-1.0 / 2.0)


m_rand = M * np.random.uniform(0.0, 1.0, N)
r_rand = r_of_m(m_rand, a, M)
phi_rand = np.random.uniform(0.0, 2 * np.pi, N)
theta_rand = np.arccos(np.random.uniform(-1.0, 1.0, N))

x = r_rand * np.sin(theta_rand) * np.cos(phi_rand)
y = r_rand * np.sin(theta_rand) * np.sin(phi_rand)
z = r_rand * np.cos(theta_rand)

X = np.array([x, y, z]).transpose()

### Generate Velocities
# Sample the DF of the Plummer model at the radii generated above (give vr,vt)
# Dejonghe (1987) Anisotropic Distribution Functionl
normalization = 3.0 * sci.gamma(6.0 - q) / (2.0 * (2.0 * np.pi) ** (2.5))
units = G ** (q - 5.0) * M ** (q - 4.0) * a ** (2.0 - q)


def H(E, L):
    x = L ** 2 / (2.0 * E * a ** 2)
    if x <= 1.0:
        return 1.0 / (sci.gamma(4.5 - q)) * sci.hyp2f1(0.5 * q, q - 3.5, 1.0, x)
    else:
        return (
            1.0
            / (sci.gamma(1.0 - 0.5 * q) * sci.gamma(4.5 - 0.5 * q) * (x ** (0.5 * q)))
            * sci.hyp2f1(0.5 * q, 0.5 * q, 4.5 - 0.5 * q, 1.0 / x)
        )


def Fq(E, L):
    if E < 0.0:
        return 0.0
    else:
        return normalization * units * E ** (3.5 - q) * H(E, L)


# Helper Functions
# Total energy: E = phi - 1/2v^2. relative potential: psi = -phi. Relative energy: -E = psi - 1/2v^2
def relative_potential(r):  # = - phi
    return G * M / np.sqrt(r ** 2 + a ** 2)


def relative_energy(r, vr, vt):
    return relative_potential(r) - 0.5 * (vr ** 2 + vt ** 2)


# N.B: Angular momentum: L = r * vt

# Convenience Function for scipy.minimize negative of Fq*vt, indep. vars passed as array
def Fq_tomin(v, r):
    return -Fq(relative_energy(r, v[0], v[1]), r * v[1]) * v[1]


# Find max of DF at given radius
def fmax(r, vmax):
    # Constraint function (defined locally, dep on r)
    def vel_constr2(v):
        return vmax ** 2 - v[0] ** 2 - v[1] ** 2

    # Initial Guess
    v0 = [0.1 * vmax, 0.2 * vmax]

    # Constraint Object
    cons = {"type": "ineq", "fun": vel_constr2}

    # Args
    args = (r,)

    # Minimize through scipy.optimize.minimize
    # fm = minimize(Fq_tomin,v0,constraints=cons,method = 'COBYLA',args=args)
    fm = minimize(
        Fq_tomin,
        v0,
        constraints=cons,
        method="SLSQP",
        args=args,
        bounds=[(0, vmax), (0, vmax)],
    )

    # Min of negative df == max of df
    return -fm.fun


# Sample vr,vt from DF at given Radius through rejection sampling
def sample_vel(r):
    # Compute max velocity (equiv. condition for E>=0)
    vmax = np.sqrt(2.0 * relative_potential(r))
    # Compute max of DF at this radius
    fm = 1.1 * fmax(r, vmax)  # 1.1 factor to be sure to include max
    while True:
        # Select candidates for vr,vt based on max full velocity
        while True:
            vr = np.random.uniform(0.0, vmax)
            vt = np.random.uniform(0.0, vmax)
            if vr ** 2 + vt ** 2 <= vmax ** 2:
                break
        # Rejection Sampling on Fq
        f = np.random.uniform(0.0, fm)
        if Fq(relative_energy(r, vr, vt), r * vt) * vt >= f:
            return vr, vt


print("Sampling velocities...")
ti = time.time()
if use_parallel:
    from multiprocessing import Pool

    with Pool(nb_threads) as p:
        vels = np.array(p.map(sample_vel, r_rand)).transpose()
else:
    from tqdm import tqdm

    vels = np.empty((2, N))
    for j, r in enumerate(tqdm(r_rand)):
        vels[:, j] = sample_vel(r)
tf = time.time()
print("Sampling took " + "{:.2f}".format(tf - ti) + " seconds.")

# Convert to Cartesian
# First: project vt on e_theta, e_phi with random orientation
alph = np.random.uniform(0, 2 * np.pi, N)
sgn = np.random.choice([-1, 1], size=N)
vphi = vels[1] * np.cos(alph)
vtheta = vels[1] * np.sin(alph)
# project vr on e_r (random sign)
vr = sgn * vels[0]

# Convert Spherical to cartesian coordinates
v_x = (
    np.sin(theta_rand) * np.cos(phi_rand) * vr
    + np.cos(theta_rand) * np.cos(phi_rand) * vtheta
    - np.sin(phi_rand) * vphi
)
v_y = (
    np.sin(theta_rand) * np.sin(phi_rand) * vr
    + np.cos(theta_rand) * np.sin(phi_rand) * vtheta
    + np.cos(phi_rand) * vphi
)
v_z = np.cos(theta_rand) * vr - np.sin(theta_rand) * vtheta

# Create velocity array
V = np.array([v_x, v_y, v_z]).transpose()

### Generate masses

m = M / N * np.ones(N)

### Exclude extreme outliers
idx = r_rand < args.bound
X = X[idx]
V = V[idx]
m = m[idx]
new_N = len(m)

### Write to hdf5
galactic_units = unyt.UnitSystem(
    "galactic",
    unyt.kpc,
    unyt.unyt_quantity(1e10, units=unyt.Solar_Mass),
    unyt.unyt_quantity(1.0, units=unyt.s * unyt.Mpc / unyt.km).to(unyt.Gyr),
)
wrt = Writer(galactic_units, args.boxsize * unyt.kpc)
wrt.dark_matter.coordinates = X * unyt.kpc
wrt.dark_matter.velocities = V * (unyt.km / unyt.s)
wrt.dark_matter.masses = m * 1e10 * unyt.msun
wrt.dark_matter.particle_ids = np.arange(new_N)
wrt.write(fname)
print("Writing IC file...")

### Optional: Pickle ic arrays for future use
if pickle_ics:
    import pickle as pkl

    with open("X.pkl", "wb") as f:
        pkl.dump(X[:], f)

    with open("V.pkl", "wb") as f:
        pkl.dump(V[:], f)

    with open("M.pkl", "wb") as f:
        pkl.dump(m[:], f)
