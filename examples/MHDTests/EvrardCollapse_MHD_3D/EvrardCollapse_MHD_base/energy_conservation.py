from swiftsimio import load_statistics
import numpy as np
import matplotlib.pyplot as plt

data = load_statistics("statistics.txt")

print(data)

Ekin = data.kin_energy.value
Eint = data.int_energy.value
Epot = data.pot_energy.value
Emag = data.mag_energy.value

plt.plot(Ekin, 'b', label='Ekin')
plt.plot(Eint, 'g', label='Eint')
plt.plot(Epot, 'y', label='Epot')
plt.plot(Emag, 'r', label='Emag')
plt.plot(Ekin+Eint+Epot, 'm', label='Ehydro')
plt.plot(Ekin+Eint+Emag+Epot, 'k', label='Etot')

plt.xlabel('Time')
plt.ylabel('Energies')
plt.xticks(np.arange(0,400, step=100), np.arange(4))
plt.legend()

plt.savefig("EnergyPlot.png", dpi=200)
