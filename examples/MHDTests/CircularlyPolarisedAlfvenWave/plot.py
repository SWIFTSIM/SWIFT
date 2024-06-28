from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt

wl = 1.0        
k = 2.0 * np.pi / wl
                                                                                                                     
v0 = 0.1
B0 = 0.1

with h5py.File(sys.argv[1], "r") as handle:
    boxsize = handle["Header"].attrs["BoxSize"][0]
    t = handle["Header"].attrs["Time"]
    
    rho = handle["PartType0/Densities"][:]

    pos = handle["PartType0/Coordinates"][:]
    v = handle["PartType0/Velocities"][:]

    h = handle["PartType0/SmoothingLengths"][:]

    B = handle["PartType0/MagneticFluxDensities"][:]
    divB = handle["PartType0/MagneticDivergences"][:]

fig, ax = plt.subplots(2,3, figsize=(11,6))

ax00b = ax[0,0].twinx()
ax10b = ax[1,0].twinx()

for filename in sys.argv[1:]:
    
    data = load(filename)
    #t = data.metadata.time
    pos = data.gas.coordinates
    x = pos[:,0] 
    v = data.gas.velocities
    B = data.gas.magnetic_flux_densities
    divB = data.gas.magnetic_divergences
    
    normB = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)
    err = np.log10(h * abs(divB) / normB) 

    for ind, axi in enumerate(ax[0,:]):
        axi.scatter(x, B[:,ind], s=0.2, label=filename)
    
    for ind, axi in enumerate(ax[1,:]):
        axi.scatter(x, v[:,ind], s=0.2, label=filename)

    ax00b.scatter(x, rho, s=0.2, c='r')   
    ax10b.scatter(x, err, s=0.2, c='r')  
  
pts = 1000

x1exact = np.linspace(0.0, max(x), pts)
vexact = [np.zeros(pts), v0 * np.sin(k * x1exact), v0 * np.cos(k * x1exact)]
Bexact = [np.ones(pts), B0 * np.sin(k * x1exact), B0 * np.cos(k * x1exact)]

for ind, axi in enumerate(ax[0,:]):
    ax[0, ind].plot(x1exact, Bexact[ind], "k-", lw=0.5, label="Exact Solution")
    ax[1, ind].plot(x1exact, vexact[ind], "k-", lw=0.5, label="Exact Solution")

#ax[0,0].set_ylim([0.99,1.01])
#ax[1,0].set_ylim([-0.01,0.01])
#ax[0,0].set_yticks([0.99, 1.00, 1.01])
#ax[1,0].set_yticks([-0.01, 0.00, 0.01])

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
