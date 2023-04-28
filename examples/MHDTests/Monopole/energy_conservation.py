from swiftsimio import load_statistics
import numpy as np
import matplotlib.pyplot as plt
import sys

#print(data)

def getEnergy(filename):
    
    data = load_statistics(filename)

    Ekin = data.kin_energy.value
    Eint = data.int_energy.value
    Emag = data.mag_energy.value

    Etot = Ekin + Eint + Emag

    Enorm = Etot[0]

    E = np.stack((Ekin, Eint, Emag, Etot), axis = 0)
    E = E/Enorm
    
    print(E.shape)

    return E

def getdivB(filename):
    
    data = load_statistics(filename)

    return data.divb_err.value
   
Econtrol = getEnergy("statistics.txt")
#EfullMHD = getEnergy("statistics_fullMHD.txt")

divBcontrol = getdivB("statistics.txt")
#divBfullMHD = getdivB("statistics_fullMHD.txt")

labels = ['Ekin', 'Eint', 'Emag', 'Etot'] 
colours = ['b', 'g', 'r', 'k']

fig, ax = plt.subplots(1, 2, figsize=(10,5))

for ii in [0,2]:

    ax[0].plot(Econtrol[ii,:], color=colours[ii], label=labels[ii])
#    ax[0].plot(EfullMHD[ii,:], color=colours[ii], linestyle=':')
    
ax[1].semilogy(divBcontrol/27648, 'c', label='Control')
#ax[1].semilogy(divBfullMHD/27648, 'c:', label='Dedner')

ax[0].set_xlabel('Time')
ax[0].set_ylabel('Energy')

#ax[0].set_ylim([-0.5, 1.5])

ax[1].set_xlabel('Time')
ax[1].set_ylabel('Divergence Error')

ax[0].grid(alpha=0.2)
ax[1].grid(alpha=0.2)

ax[0].legend()
ax[1].legend()

plt.tight_layout()

plt.savefig(sys.argv[1], dpi=200)
