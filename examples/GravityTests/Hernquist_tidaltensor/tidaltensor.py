import swiftsimio
import matplotlib.pyplot as plt
import unyt as u
import numpy as np
import sys
# import re
# import glob

G = swiftsimio.cosmo_quantity(1, u.G, comoving=True, scale_factor=1, scale_exponent=3)

def rho(r, M, a, mres):
    '''
    Preffered units:
    r in kpc
    M in Msun
    a in kpc
    mres in Msun
    '''
    return  M / (2 * np.pi) * a / (r * (r + a) ** 3)

def shell_vol(rin, rout):
    ''' Volume of spherical shell between radii rin and rout'''

    return 4 * np.pi / 3 * (rout ** 3 - rin ** 3)

def ten_norm_func(r, M, a):
    return G * M * np.sqrt(6 * r ** 2 + 4 * r * a + 2 * a ** 2) / (r * (r + a) ** 3)

def running_avg(shell_mask, arr):
    N_shells = len(shell_mask)
    masked = shell_mask * np.array([arr,] * N_shells)
    masked[masked == 0] = np.nan
    return np.nanmean(masked, axis=1)

snapfile = 'output_0000.hdf5'
snap_mask = swiftsimio.mask(snapfile)
snap = swiftsimio.load(snapfile)
G = swiftsimio.cosmo_quantity(1, u.G, comoving=True, scale_factor=1, scale_exponent=3)

# Read in some data
x_gas = snap.gas.coordinates
x_cen = np.median(x_gas, axis=0)
x_gas -= x_cen                      # center the data around median
r_gas = np.sqrt(np.sum(x_gas * x_gas, axis=1)).to(u.kpc)
dens = snap.gas.densities.to(u.mh/u.cm**3)
epsilon_sim = snap.gas.softenings.to(u.kpc)

# Calculate the normed value of the tensor from data
ten2 = snap.gas.tidal_tensors ** 2
ten_norm = np.sqrt(np.sum(ten2, axis=1) + np.sum(ten2[:, 3:], axis=1)).to(u.Myr ** -2)

# Hernquist halo characteristics (from IC file)
M = 1e10 * u.msun # Msun
a = 10  * u.kpc # kpc
mres = 1e5 * u.msun

# Calculate expected tidal tensor and softening based on IC file:
r_ = np.logspace(-1, 2.5, num=50) * u.kpc
ten_ = ten_norm_func(r_, M, a).to(u.Myr ** -2)
epsilon_ = ((G * mres / ten_) ** (1/3)).to(u.kpc)


# Calculate running averages:
shell_mask = np.zeros([49, len(r_gas)])
for i, (rmin, rmax) in enumerate(zip(r_[:-1], r_[1:])):
    shell_mask[i] = (r_gas > rmin) & (r_gas < rmax)

avg_dens = running_avg(shell_mask, dens)
avg_ten = running_avg(shell_mask, ten_norm)
avg_epsilon = running_avg(shell_mask, epsilon_sim)



# Plot..
fig, axs = plt.subplots(3, 1, figsize=[5, 12])
fig.subplots_adjust(hspace=0.3)
fig.suptitle(r'Hernquist halo M = $10^{11}$M$_\odot$, c = 10 kpc' + \
    '\n0th snap and theoretical profiles', y=0.95)

ax = axs[0]
ax.scatter(r_gas, dens, zorder=1, alpha=0.2, s=1)
ax.plot(r_, rho(r_, M, a, mres).to(u.mh/u.cm**3), c='k', label='Hernquist profile')
ax.plot((r_[1:] + r_[:-1]) / 2, avg_dens, linestyle='--', c='r', label='Running average')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Distance [kpc]')
ax.set_ylabel(r'$\rho$ [$m_h/$cm$^3$]')
ax.set_title('Particle densities')
ax.set_xlim(1e-1, 1e2)
legend = ax.legend(loc='lower left')

ax = axs[1]
ax.scatter(r_gas, ten_norm, alpha=0.3, label='sim', s=1)
ax.plot((r_[1:] + r_[:-1]) / 2, avg_ten, linestyle='--', c='r', label='Running average')
ax.plot(r_, ten_, c='k')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Distance [kpc]')
ax.set_ylabel(r'|T| [Myr$^{-2}$]')
ax.set_title('Particle Tidal tensor')
ax.set_xlim(1e-1, 1e2)

ax = axs[2]
ax.scatter(r_gas, epsilon_sim, alpha=0.3, label='sim', s=1)
ax.plot(r_, epsilon_, c='k')
ax.plot((r_[1:] + r_[:-1]) / 2, avg_epsilon, linestyle='--', c='r')
ax.axhline(0.2, c='dimgrey', alpha=0.8, label='Default softening')
# ax.axhline(0.01, c='dimgrey', alpha=0.8, linestyle='--')
ax.set_xlabel('Distance [kpc]')
ax.set_ylabel(r'$\epsilon$ [kpc]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Softening length')
ax.set_xlim(1e-1, 1e2)


n = 'particleTT.png'
plt.savefig(f'{n}', bbox_inches='tight')
print(f'Saving: {n}')

