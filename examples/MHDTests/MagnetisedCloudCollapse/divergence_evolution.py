
#from swiftsimio import load_statistics
from unyt import g, s, statA
import numpy as np
import matplotlib.pyplot as plt
import sys

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

file = sys.argv[1]

data = np.loadtxt(file)

t = data[:,1]
   
divB_err = data[:,35] 
Brms = data[:,38]

ln1 = ax1.semilogy(t, Brms/Brms[0],'k-', label=r'$B_{rms}$')
ln2 = ax2.semilogy(t, divB_err, 'r:', label=r'Error')
    
ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r'$B_{rms} / B_{rms} (t=0)$')
ax2.set_ylabel(r'$Error$')

ax1.set_xlim([0.0, 0.05])

lns = ln1 + ln2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs)

plt.tight_layout()

plt.savefig('DivergencePlot.png', dpi=200)
