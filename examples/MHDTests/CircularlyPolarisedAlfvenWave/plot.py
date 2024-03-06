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

sinb = 0.0 # 2.0 / np.sqrt(5) 
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
filename = sys.argv[1]

# Read in snapshot (temporary, once we like the IC file we'll go back to using swiftsimio)
with h5py.File(filename, "r") as handle:
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"]
    
    rho = handle["PartType0/Densities"][:]

    pos = handle["PartType0/Coordinates"][:]
    v = handle["PartType0/Velocities"][:]

    h = handle["PartType0/SmoothingLengths"][:]

    B = handle["PartType0/MagneticFluxDensities"][:]
    divB = handle["PartType0/MagneticDivergences"][:]

# Plot things
fig, ax = plt.subplots(2,3, figsize=(11,6))

ax00b = ax[0,0].twinx()
ax10b = ax[1,0].twinx()
 
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

for ind, axi in enumerate(ax[0,:]):
    axi.scatter(posp[0,:], Bp[ind,:], s=0.2, label=filename)
  
for ind, axi in enumerate(ax[1,:]):
    axi.scatter(posp[0,:], vp[ind,:], s=0.2, label=filename)

ax00b.scatter(posp[0,:], rho, s=0.2, c='g')   
ax10b.scatter(posp[0,:], err, s=0.2, c='r')  
  
pts = 1000

x1exact = np.linspace(0.0, max(posp[0,:]), pts)
vexact = [np.zeros(pts), v0 * np.sin(k * x1exact), v0 * np.cos(k * x1exact)]
Bexact = [np.ones(pts), B0 * np.sin(k * x1exact), B0 * np.cos(k * x1exact)]

for ind, axi in enumerate(ax[0,:]):
    ax[0, ind].plot(x1exact, Bexact[ind], "k-", lw=0.5, label="Exact Solution")
    ax[1, ind].plot(x1exact, vexact[ind], "k-", lw=0.5, label="Exact Solution")

labelsv = [r'$v_1$', r'$v_2$', r'$v_3$']
labelsB = [r'$B_1$', r'$B_2$', r'$B_3$']

for ind, axi in enumerate(ax[0,:]):
    ax[1,ind].set_xlabel(r'$x_1$')
    ax[0,ind].set_ylabel(labelsB[ind])
    ax[1,ind].set_ylabel(labelsv[ind])
    ax[0,ind].set_xticks([])

ax00b.set_ylabel(r'$\rho$')
ax10b.set_ylabel(r'$err_{\nabla \cdot B}$')

plt.tight_layout()

plt.savefig("AlfvenWaves.png", dpi=100)
