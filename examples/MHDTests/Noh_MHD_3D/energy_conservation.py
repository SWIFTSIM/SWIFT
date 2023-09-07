from swiftsimio import load_statistics
import numpy as np
import matplotlib.pyplot as plt

data = load_statistics("statistics.txt")

print(data)

Ekin = data.kin_energy.value
Eint = data.int_energy.value
Epot = data.pot_energy.value
Emag = data.mag_energy.value

divErr = data.divb_err.value

plt.plot(Ekin, "b", label="Ekin")
plt.plot(Eint, "g", label="Eint")
plt.plot(Emag, "r", label="Emag")
# plt.plot(Ekin+Eint, 'm', label='Ehydro')
plt.plot(Ekin + Eint + Emag, "k", label="Etot")

plt.plot(divErr, "c", label="Error")

plt.xlabel("Time")
plt.ylabel("Energies")
# plt.xticks(np.arange(0,400, step=100), np.arange(4))
# plt.xlim([0,420])
plt.ylim([-0.5, 4.5])
plt.legend()

plt.savefig("EnergyPlot.png", dpi=200)
