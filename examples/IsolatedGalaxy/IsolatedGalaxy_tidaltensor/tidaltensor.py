import swiftsimio
import matplotlib.pyplot as plt
import unyt as u
import numpy as np
import sys
import re
import glob

G = swiftsimio.cosmo_quantity(1, u.G, comoving=True, scale_factor=1, scale_exponent=3)

def ten_norm_func(r, M, a):
    ''' Theoretical tidal tensor norm for a Hernquist potential.'''

    return G * M * np.sqrt(6 * r ** 2 + 4 * r * a + 2 * a ** 2) / (r * (r + a) ** 3)


snapfile = './output_0000.hdf5'

snap = swiftsimio.load(snapfile)

# Read in some data
x_dm = snap.dark_matter.coordinates
x_cen = np.median(x_dm, axis=0)
x_dm -= x_cen                      # center the data around median
r_dm = np.sqrt(np.sum(x_dm * x_dm, axis=1)).to(u.kpc)

x_stars = snap.stars.coordinates
x_cen = np.median(x_stars, axis=0)
x_stars -= x_cen                      # center the data around median
r_stars = np.sqrt(np.sum(x_stars * x_stars, axis=1)).to(u.kpc)

# Calculate the normed value of the tensor from data
ten2 = snap.dark_matter.tidal_tensors ** 2
ten_norm_dm = np.sqrt(np.sum(ten2, axis=1) + np.sum(ten2[:, 3:], axis=1)).to(u.Myr ** -2)
ten2 = snap.stars.tidal_tensors ** 2
ten_norm_stars = np.sqrt(np.sum(ten2, axis=1) + np.sum(ten2[:, 3:], axis=1)).to(u.Myr ** -2)

# Halo parameters
M = 1.99e12 * u.msun
a = 52.19 * u.kpc
hernquist_soft = 0.2 * u.kpc # gravitational softening


# Calculate the SOFTENEND tidal tensor profile
r_ = np.logspace(-1, 2, num=15) * u.kpc
ten_ = ten_norm_func(np.sqrt(r_ ** 2  + hernquist_soft ** 2), M, a).to(u.Myr ** -2)

# and the unsoftened one just for fun
ten_unsoft = ten_norm_func(r_, M, a).to(u.Myr ** -2)

# Plot...
fig, ax = plt.subplots(figsize=[5,4])
plt.tight_layout()
ax.scatter(r_dm, ten_norm_dm, c='C1', alpha=0.5, label='DM', s=1)
ax.scatter(r_stars, ten_norm_stars, c='C0', alpha=0.5, label='stars', s=1)
ax.plot(r_, ten_, c='k', linewidth=1, label='softened profile')
ax.plot(r_, ten_unsoft, c='k', linestyle=':', alpha=0.5, label='unsoftened profile')
ax.axvline(0.2, label=r'Particle softening $\varepsilon$', c='k', linestyle='--')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Distance [kpc]')
ax.set_ylabel(r'|T| [Myr$^{-2}$]')
ax.set_title('Particle Tidal tensor')
ax.set_xlim(1e-1, 1e2)
ax.set_ylim(bottom=1e-5)
ax.legend()


n = 'particleTT.png'
plt.savefig(f'{n}', bbox_inches='tight')
print(f'Saving: {n}')


