from swiftsimio import load
import numpy as np
import sys
import matplotlib.pyplot as plt

fig = plt.figure()

gs = fig.add_gridspec(5, hspace=0)
axs = gs.subplots(sharex=True)

colours = ['b', 'r', 'g', 'm', 'c']

x1exact = np.linspace(0.0, 3.0, 1000)
rhoexact = np.ones(1000)
B1exact = np.ones(1000)
B2exact = 0.1 * np.sin(2 * np.pi * x1exact)
B3exact = 0.1 * np.cos(2 * np.pi * x1exact)

axs[0].plot(x1exact, rhoexact, 'k--', lw=0.5, label='Exact Solution')
axs[1].plot(x1exact, B1exact, 'k--', lw=0.5)
axs[2].plot(x1exact, B2exact, 'k--', lw=0.5)
axs[3].plot(x1exact, B3exact, 'k--', lw=0.5)

for colour, file in zip(colours, sys.argv[1:]):

    # Get particle attribuets
    data = load(file)

    print(data)

    pos = data.gas.coordinates
    h = data.gas.smoothing_lengths
    rho = data.gas.densities
    B = data.gas.magnetic_flux_densities
    normB = np.sqrt(B[:,0] * B[:,0] + B[:,1] * B[:,1] + B[:,2] * B[:,2]) 
    divB = data.gas.magnetic_divergences
    
    x1 = (pos[:, 0] + 2 * pos[:, 1] + 2 * pos[:, 2]) / 3
    B1 = (B[:, 0] + 2 * B[:, 1] + 2 * B[:, 2]) / 3
    B2 = (B[:, 1] - 2 * B[:, 0]) / np.sqrt(5)
    B3 = (- 2 * B[:, 0] - 4 * B[:, 1] + 5 * B[:, 2]) / (3 * np.sqrt(5)) 

    sort_ind = np.argsort(x1)

    err = h[sort_ind] * abs(divB[sort_ind]) / normB[sort_ind]
    
    axs[0].plot(x1[sort_ind], rho[sort_ind], c=colour, lw=0.5, label=file)
    axs[1].plot(x1[sort_ind], B1[sort_ind], c=colour, lw=0.5, label=file)
    axs[2].plot(x1[sort_ind], B2[sort_ind], c=colour, lw=0.5, label=file)
    axs[3].plot(x1[sort_ind], B3[sort_ind], c=colour, lw=0.5, label=file)
    axs[4].semilogy(x1[sort_ind], err, c=colour, ls=':', lw=0.5, label=file)    

for ax in axs:
    ax.label_outer()

axs[3].set_xlabel(r"$x_1$")
axs[0].set_ylabel(r"$\rho$")
axs[1].set_ylabel(r"$B_1$")
axs[2].set_ylabel(r"$B_2$")
axs[3].set_ylabel(r"$B_3$")
axs[4].set_ylabel(r"$h \: \nabla \cdot B / |B|$")

axs[4].legend()

plt.savefig("CircularlyPolarisedAlfvenWave.png", dpi=200)
