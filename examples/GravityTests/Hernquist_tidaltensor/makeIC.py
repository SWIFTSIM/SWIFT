# Adapted from swiftsimio 'Creating Initial Conditions' tutorial


import numpy as np
import matplotlib.pyplot as plt
import swiftsimio
import unyt as u


def rho(r, M, a, mres):
    '''
    Preferred units:
    r in kpc
    M in Msun
    a in kpc
    mres in Msun
    '''
    return  M * a / (2 * np.pi * r * (r + a) ** 3)

def shell_vol(rin, rout):
    return 4 * np.pi / 3 * (rout ** 3 - rin ** 3)

def convert_to_xyz(rad, angles):
    # angles = N x 2: [N, theta/phi]
    loc = np.zeros((len(rad), 3))

    loc[:, 0] = rad * np.sin(angles[:, 0]) * np.sin(angles[:, 1])
    loc[:, 1] = rad * np.sin(angles[:, 0]) * np.cos(angles[:, 1])
    loc[:, 2] = rad * np.cos(angles[:, 0])

    return loc

a = 10 # kpc
M = 1e10 # Msun
mres = 1e5 # Msun


shell_edges = np.logspace(-1, 3)
shell_width = (shell_edges[1:] + shell_edges[:-1]) / 2
shell_cen = np.sqrt(shell_edges[:-1] * shell_edges[1:])

rho_edges = rho(shell_edges, M, a, mres)
rho_cen = (rho_edges[1:] + rho_edges[:-1]) / 2

loc_all = []

# Distribute particles randomly in each shell.
print('Seed positions...')
for i, r in enumerate(shell_cen):
    V = shell_vol(shell_edges[i], shell_edges[i+1])
    # No. particles needed to have desired density
    N = int(rho_cen[i] * V / mres)

    # seed theta, phi
    angles = np.random.uniform(0, 1, (N, 2))
    angles[:, 1] *= (2 * np.pi)
    angles[:, 0] = np.arccos(2 * angles[:, 0] - 1)

    rad = np.random.uniform(shell_edges[i], shell_edges[i+1], N)
    loc = convert_to_xyz(rad, angles)
    loc += 2000         # offset postitions 

    loc_all.append(loc)
loc_all = np.vstack(loc_all)

n_p = len(loc_all)
print('Total particles:' , n_p)

# WRITE
# ------
print('Creating IC file...')

# Boxsize is 4 x 4 x 4 Mpc
lbox = 4000
boxsize = swiftsimio.cosmo_array(
     [lbox, lbox, lbox],
     u.kpc,
     comoving=True,
     scale_factor=1.0,
     scale_exponent=1,
)


# Galactic unit system: kpc, Msun, kyr
fw = swiftsimio.Writer(unit_system='galactic', boxsize=boxsize)
# fw.gas.coordinates = loc_all * u.kpc

fw.dark_matter.coordinates = swiftsimio.cosmo_array(
    loc_all,
    u.kpc,
    comoving=True,
    scale_factor=fw.scale_factor,
    scale_exponent=1,
)

# Random velocities from 0 to 1 km/s
fw.dark_matter.velocities = swiftsimio.cosmo_array(
    np.random.rand(n_p, 3),
    u.km / u.s,
    comoving=True,
    scale_factor=fw.scale_factor,
    scale_exponent=1,
)

# Generate uniform masses as 10^5 solar masses for each particle
fw.dark_matter.masses = swiftsimio.cosmo_array(
    np.ones(n_p, dtype=float) * 1e5,
    u.msun,
    comoving=True,
    scale_factor=fw.scale_factor,
    scale_exponent=0,
)

# # Generate internal energy corresponding to 10^4 K
# fw.gas.internal_energy = swiftsimio.cosmo_array(
#     np.ones(n_p, dtype=float) * 1e4 / 1e5,
#     u.kb * u.K / u.msun,
#     comoving=True,
#     scale_factor=fw.scale_factor,
#     scale_exponent=-2,
# )

# # fw.gas.generate_smoothing_lengths()


# # # Generate initial guess for smoothing lengths based on MIPS
# # # fw.gas.smoothing_length = (np.ones(n_p) * u.kpc)
# fw.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
# # # fw.gas.smoothing_length = 519

# print(fw.gas.smoothing_length)

# If IDs are not present, this automatically generates
fw.write("HernquistIC.hdf5")
