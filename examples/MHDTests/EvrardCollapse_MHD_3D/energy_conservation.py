from swiftsimio import load_statistics
import numpy as np
import matplotlib.pyplot as plt

data = load_statistics("statistics.txt")

print(data)

Ekin = data.kin_energy.value
Eint = data.int_energy.value
Epot = data.pot_energy.value
Emag = data.mag_energy.value

plt.semilogy(Ekin, 'b', label='Ekin')
plt.semilogy(Eint, 'g', label='Eint')
plt.semilogy(Epot, 'y', label='Epot')
plt.semilogy(Emag, 'r', label='Emag')
plt.semilogy(Ekin+Eint+Epot, 'm', label='Ehydro')
plt.semilogy(Ekin+Eint+Emag+Epot, 'k', label='Etot')

plt.xlabel('Time')
plt.ylabel('Energies')
#xplt.ylim([-2.5, 2.5])
#plt.xticks(np.arange(0,400, step=100), np.arange(4))
plt.legend()

plt.savefig("EnergyPlot.png", dpi=200)
