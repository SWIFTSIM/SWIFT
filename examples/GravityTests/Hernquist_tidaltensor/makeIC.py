
import numpy as np
import matplotlib.pyplot as plt
# from astropy import units
# import h5py as h5
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
shell_cen = shell_edges[:-1] + shell_width

rho_edges = rho(shell_edges, M, a, mres)
rho_cen = (rho_edges[1:] + rho_edges[:-1]) / 2

loc_all = []

print('Seed positions...')
for i, r in enumerate(shell_cen):
    V = shell_vol(shell_edges[i], shell_edges[i+1])
    N = int(rho_cen[i] * V / mres)

    # seed theta, phi
    angles = np.random.uniform(0, 1, (N, 2))
    angles[:, 1] *= (2 * np.pi)
    angles[:, 0] = np.arccos(2 * angles[:, 0] - 1)

    # angles = np.random.uniform(0, np.pi, (N, 2))
    # angles[:, 1] *= 2

    rad = np.random.uniform(shell_edges[i], shell_edges[i+1], N)
    loc = convert_to_xyz(rad, angles)
    loc += 2000         # offset postitions 
    # print(r, N)

    # if i == 20:
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')
    #     ax.scatter3D(loc[:, 0], loc[:, 1], loc[:, 2], s=2)
    #     plt.savefig('ball.png')
    #     exit()

    loc_all.append(loc)
loc_all = np.vstack(loc_all)

n_p = len(loc_all)
print('tot gas part' , n_p)

# PLOT TEST
# ---------
# r_ = np.linspace(1, rmax)
# r_ = np.logspace(0, 2.5)

# r_n = np.sqrt(np.sum(loc_all ** 2, axis=1))
# # print(r_n)

# hist, bin_edges = np.histogram(np.log10(r_n), bins=20)
# bin_edges = 10 ** bin_edges
# cen = bin_edges[:-1] + (bin_edges[1] - bin_edges[0]) / 2
# rho_n = hist * mres / shell_vol(bin_edges[:-1], bin_edges[1:])

# rho_ = rho(cen, M, a, mres)
# # print(rho_ / rho_n)

# # ax.hist(r_gas, bins=200)

# plt.scatter(cen, rho_n, label='Resulting density')
# plt.plot(cen, rho_, c='orange', label='Hernquist profile')
# plt.title('Expected and resulting density')
# plt.xlabel('Distance [kpc]')
# plt.ylabel(r'Density [$M_\odot/$kpc$^3$]')
# plt.xscale('log')
# plt.yscale('log')

# plt.savefig(f'Hernquist_init.png')

# exit()


# WRITE
# ------
print('Creating IC file...')

# Box is 400 kpc Mpc
boxsize = 4000 * u.kpc


# Galactic unit system: kpc, Msun, kyr
fw = swiftsimio.Writer(unit_system='galactic', boxsize=boxsize)



fw.gas.coordinates = loc_all * u.kpc



# Random velocities from 0 to 1 km/s
fw.gas.velocities = np.random.rand(n_p, 3) * (u.km / u.s)

# Generate uniform masses as 10^5 solar masses for each particle
fw.gas.masses = np.ones(n_p, dtype=float) * (1e5 * u.msun)

# Generate internal energy corresponding to 10^4 K
fw.gas.internal_energy = (
    np.ones(n_p, dtype=float) * (1e4 * u.kb * u.K) / (1e5 * u.msun)
)

# Generate initial guess for smoothing lengths based on MIPS
# fw.gas.smoothing_length = (np.ones(n_p) * u.kpc)
fw.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
fw.gas.smoothing_length /=519

# print(fw.gas.smoothing_length)

# If IDs are not present, this automatically generates
fw.write("HernquistIC.hdf5")
