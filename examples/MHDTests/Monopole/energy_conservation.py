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

    #Enorm = Etot[0]

    E = np.stack((Ekin, Eint, Emag, Etot), axis = 0)
    #E = E/np.expand_dims(E[:,0],axis=1)
    
    print(E.shape)

    return E

def getdivB(filename):
    
    data = load_statistics(filename)

    return data.divb_err.value
   
Econtrol    = getEnergy("statistics_no_cleaning.txt")
Ehyperbolic = getEnergy("statistics_hyperbolic.txt") 
EfullMHD    = getEnergy("statistics_fullMHD.txt")

divBcontrol    = getdivB("statistics_no_cleaning.txt")
divBhyperbolic = getdivB("statistics_hyperbolic.txt") 
divBfullMHD    = getdivB("statistics_fullMHD.txt")

labels = ['Ekin', 'Eint', 'Emag', 'Etot'] 
colours = ['b', 'g', 'r', 'k']

fig, ax = plt.subplots(2, 3, sharex=True, figsize=(12,6))

for ii in range(2):
    for jj in range(2):
        idx = 2*ii+jj
        ax[ii,jj].semilogy(Econtrol[idx,:], color=colours[idx]) #, label=labels[idx])
        ax[ii,jj].semilogy(Ehyperbolic[idx,:], color=colours[idx], linestyle='--')
        ax[ii,jj].plot(EfullMHD[idx,:], color=colours[idx], linestyle=':')
        ax[ii,jj].set_ylabel(labels[idx])
        #ax[ii,jj].legend()

ax[1,0].axhline(y=1/(2*np.pi), color='k', linestyle='--')
#ax[0].axhline(y=0.331008, color='k', linestyle='-')
    
ax[0,2].semilogy(divBcontrol/2304, 'c', label='Control')
ax[0,2].semilogy(divBhyperbolic/2304, 'c--', label='Hyperbolic')
ax[0,2].semilogy(divBfullMHD/2304, 'c:', label='Dedner')

#ax[0].set_xlabel('Time')
#ax[0].set_ylabel('Energy')

#ax[0].set_ylim([-0.5, 1.5])

#ax[1].set_xlabel('Time')
ax[0,2].set_ylabel('Divergence Error')

#ax[0].grid(alpha=0.2)
#ax[1].grid(alpha=0.2)

#ax[0].legend()
ax[0,2].legend()

ax[1,2].axis('off')

plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.tight_layout()

plt.savefig(sys.argv[1], dpi=200)
