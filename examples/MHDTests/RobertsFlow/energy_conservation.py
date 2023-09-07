from swiftsimio import load_statistics
import numpy as np
import matplotlib.pyplot as plt
import sys

def runDict(filename):

    data = load_statistics(filename)
                                                                                   
    Bms  = data.bms.value
    Brms = np.sqrt(Bms)
    Brms = Brms / Brms[0]

    Px = data.mom_x.value
    Py = data.mom_y.value
    Pz = data.mom_z.value

    Lx = data.ang_mom_x.value
    Ly = data.ang_mom_y.value
    Lz = data.ang_mom_z.value
    
    divB = abs(data.divb_err.value)
    P    = np.sqrt(Px**2 + Py**2 +Pz**2)
    L    = np.sqrt(Lx**2 + Ly**2 +Lz**2)

    return {"Brms":Brms, "divB":divB, "Pressure":P, "AngMom":L}

fig, ax = plt.subplots(2, 2, figsize=(8,8))

for file in sys.argv[1:]:
    
    print(file)

    stats = runDict(file)
    
    ax[0,0].semilogy(stats["Brms"], label=file)
    ax[0,1].semilogy(stats["divB"]/65536)
    ax[1,0].semilogy(stats["Pressure"])
    ax[1,1].semilogy(stats["AngMom"])

t_anal = np.linspace(0,25,1000)
Brms_anal = np.exp(0.05*(t_anal + 25))

ax[0,0].semilogy(t_anal, Brms_anal, 'k--')

ax[0,0].set_ylabel(r'$B_{<rms>} / B_{<rms>} (t=0)$')
ax[0,0].legend(loc="best")
    
ax[0,1].set_ylabel('Divergence Error')
ax[0,1].set_ylim([1e-4,10])
ax[1,0].set_ylabel('Linear momentum')

ax[1,1].set_ylabel('Angular momentum')

'''
tff  = 3e4 * 3.156e7
tsim = 0.05 * 3.08567758149E13

locs   = [0, 100, 200, 300, 400, 500]
labels = tsim / tff * np.linspace(0.0, 1.0, 6)
labels = np.round(labels, decimals=2)

for axi in ax:
    for axii in axi:
        axii.set_xlabel(r'Time $[t_{ff}]$')
        axii.set_xticks(locs, labels)
'''

plt.tight_layout()

plt.savefig("EnergyPlot", dpi=200)
