from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

# Test parameters 
wl = 1.0        
k = 2.0 * np.pi / wl
                                                                                                                     
v0 = 0.1
B0 = 0.1

sina = 0.0 # 2.0 / 3.0 
cosa = np.sqrt(1 - sina ** 2)

sinb = 0.5 # 2.0 / np.sqrt(5) 
cosb = np.sqrt(1 - sinb ** 2)

# Rotation matrix to rotate frame S to the projected direction frame Sp ...
Rotation_Sp_to_S = np.array(
    [
        [cosa * cosb, -sinb, -sina * cosb],
        [cosa * sinb, cosb, -sina * sinb],
        [sina, 0, cosa],
    ]
)

# ... and its inverse to go from Sp to S
Rotation_S_to_Sp = Rotation_Sp_to_S.T

# File to read
snap_num = int(sys.argv[1])
filename = "CircularlyPolarisedAlfvenWave_%0.4d.hdf5"%snap_num

# Read in snapshot (temporary, once we like the IC file we'll go back to using swiftsimio)
with h5py.File(filename, "r") as handle:
    boxsize = handle["Header"].attrs["BoxSize"]
    t = handle["Header"].attrs["Time"]

    rho = handle["PartType0/Densities"][:]
    P = handle["PartType0/Pressures"][:]
    u = handle["PartType0/InternalEnergies"][:]
    pos = handle["PartType0/Coordinates"][:]
    ids = handle["PartType0/ParticleIDs"][:]
    v = handle["PartType0/Velocities"][:]
    h = handle["PartType0/SmoothingLengths"][:]
    B = handle["PartType0/MagneticFluxDensities"][:]
    divB = handle["PartType0/MagneticDivergences"][:]

print(np.min(u), np.max(u))
    
# Plot things
fig, ax = plt.subplots(2, 4, figsize=(13,6))

fig.suptitle('t=%.3f'%t, fontsize=16)

#ax00b = ax[0,0].twinx()
#ax10b = ax[1,0].twinx()
 
'''
data = load(filename)
t = data.metadata.time
pos = data.gas.coordinates.T
v = data.gas.velocities.T
B = data.gas.magnetic_flux_densities.T
'''
    
normB = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)
err = np.log10(h * abs(divB) / normB) 

pos = pos.T
v = v.T
B = B.T

posp = np.dot(Rotation_S_to_Sp, pos)
vp = np.dot(Rotation_S_to_Sp, v)
Bp = np.dot(Rotation_S_to_Sp, B)

for ind, axi in enumerate(ax[0,:3]):
    axi.scatter(posp[0,:], Bp[ind,:], s=0.2, label=filename)
  
for ind, axi in enumerate(ax[1,:3]):
    axi.scatter(posp[0,:], vp[ind,:], s=0.2, label=filename)

ax[0,3].scatter(posp[0,:], rho, s=0.2)
ax[1,3].scatter(posp[0,:], P, s=0.2)

# Exact solution
pts = 1000
x1exact = np.linspace(0.0, max(posp[0,:]), pts)
vexact = [np.zeros(pts), v0 * np.sin(k * x1exact), v0 * np.cos(k * x1exact)]
Bexact = [np.ones(pts), B0 * np.sin(k * x1exact), B0 * np.cos(k * x1exact)]
Pexact = np.ones(pts) * 0.1
rhoexact = np.ones(pts) * 1

for ind, axi in enumerate(ax[0,:3]):
    ax[0, ind].plot(x1exact, Bexact[ind], "k-", lw=0.5, label="Exact Solution")
    ax[1, ind].plot(x1exact, vexact[ind], "k-", lw=0.5, label="Exact Solution")

ax[0,3].plot(x1exact, rhoexact, "k-", lw=0.5, label="Exact Solution")
ax[1,3].plot(x1exact, Pexact, "k-", lw=0.5, label="Exact Solution")
    
labelsv = [r'$v_{par}$', r'$v_{per}$', r'$v_z$']
labelsB = [r'$B_{par}$', r'$B_{per}$', r'$B_z$']

for ind in range(3):
    ax[1,ind].set_xlabel(r'$x_{par}$')
    ax[0,ind].set_ylabel(labelsB[ind])
    ax[1,ind].set_ylabel(labelsv[ind])
    ax[0,ind].set_xticks([])
    ax[0,ind].set_ylim(-0.11,0.11)
    ax[1,ind].set_ylim(-0.11,0.11)
    ax[0,ind].set_xlim(0,2)
    ax[1,ind].set_xlim(0,2)

ax[1,3].set_xlabel(r'$x_{par}$')
ax[0,3].set_ylabel(r'$\rho$')
ax[0,3].set_ylim(0.9,1.1)

ax[1,3].set_xlabel(r'$x_{par}$')
ax[1,3].set_ylabel(r'$P$')
ax[1,3].set_ylim(0.,0.2)
    
ax[0,0].set_ylim(0.9,1.1)

ax[0, 1].hlines(0.1, -5, 5, ls=':', lw=0.5, color='k')
ax[0, 1].hlines(-0.1, -5, 5, ls=':', lw=0.5, color='k')
ax[1, 1].hlines(0.1, -5, 5, ls=':', lw=0.5, color='k')
ax[1, 1].hlines(-0.1, -5, 5, ls=':', lw=0.5, color='k')


plt.tight_layout()

plt.savefig("AlfvenWaves_%0.4d.png"%snap_num, dpi=100)

################################################################33

my_id = 473

# Plot the particle positions

fig, ax = plt.subplots(1,1, figsize=(5,8))
fig.suptitle('t=%.3f'%t, fontsize=16)

mask = ids == my_id

ax.scatter(pos.T[:,0], pos.T[:,1], s=3)
ax.scatter(pos.T[mask,0], pos.T[mask,1], s=30, color='r')
ax.hlines(0, -5, 5, ls='--', lw=0.5, color='k')
ax.hlines(boxsize[1], -5, 5, ls='--', lw=0.5, color='k')
ax.vlines(0, -5, 5, ls='--', lw=0.5, color='k')
ax.vlines(boxsize[0], -5, 5, ls='--', lw=0.5, color='k')

ax.set_xlim(-0.1, boxsize[0]+0.1)
ax.set_ylim(-0.1, 2.1)

plt.tight_layout()

plt.savefig("Positions_%0.4d.png"%snap_num, dpi=100)
