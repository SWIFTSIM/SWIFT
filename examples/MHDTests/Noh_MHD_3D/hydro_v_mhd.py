from swiftsimio import load_statistics
import numpy as np
import matplotlib.pyplot as plt

data_hydro = load_statistics("statistics_hydro.txt")
data_mhd = load_statistics("statistics_mhd.txt")

Ekin_hydro = data_hydro.kin_energy.value
Eint_hydro = data_hydro.int_energy.value
Epot_hydro = data_hydro.pot_energy.value
Emag_hydro = data_hydro.mag_energy.value

Ekin_mhd = data_mhd.kin_energy.value
Eint_mhd = data_mhd.int_energy.value
Epot_mhd = data_mhd.pot_energy.value
Emag_mhd = data_mhd.mag_energy.value

plt.plot((Ekin_mhd - Ekin_hydro[: int(len(Ekin_mhd))]) / Ekin_mhd, "b", label="Ekin")
plt.plot((Eint_mhd - Eint_hydro[: int(len(Ekin_mhd))]) / Eint_mhd, "g", label="Eint")
plt.plot((Emag_mhd - Emag_hydro[: int(len(Ekin_mhd))]) / Emag_mhd, "r", label="Emag")
# plt.plot(Ekin+Eint, 'm', label='Ehydro')
# plt.plot(Ekin+Eint+Emag, 'k', label='Etot')

plt.xlabel("Time")
plt.ylabel("MHD - hydro relative energy diff")
# plt.xticks(np.arange(0,400, step=100), np.arange(4))
plt.xlim([0, 420])
plt.legend()

plt.savefig("EnergyDiff.png", dpi=200)
