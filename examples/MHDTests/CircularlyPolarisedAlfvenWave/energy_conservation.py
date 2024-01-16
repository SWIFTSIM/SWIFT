from swiftsimio import load_statistics
import matplotlib.pyplot as plt

data = load_statistics("statistics.txt")

Ekin = data.kin_energy.value
Eint = data.int_energy.value
Emag = data.mag_energy.value

divB_err = data.divb_err.value

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

ax1.plot(Ekin, "b", label="Ekin")
ax1.plot(Eint, "g", label="Eint")
ax1.plot(Emag, "r", label="Emag")
ax1.plot(Ekin + Eint + Emag, "k", label="Etot")
ax2.plot(divB_err, "c", label="Error")

ax1.set_xlabel("Time")
ax1.set_ylabel("Energy")
ax2.set_ylabel("Divergence Error")

ax1.legend()

plt.savefig("EnergyPlot.png", dpi=200)
