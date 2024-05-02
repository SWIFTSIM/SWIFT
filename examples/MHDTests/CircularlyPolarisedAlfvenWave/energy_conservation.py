from swiftsimio import load_statistics
import matplotlib.pyplot as plt

data = load_statistics("statistics.txt")

time = data.time
Ekin = data.kin_energy.value
Eint = data.int_energy.value
Emag = data.mag_energy.value

divB_err = data.divb_err.value

fig, ax1 = plt.subplots()

#ax2 = ax1.twinx()

Etot = Ekin + Eint + Emag

ax1.plot(time, Ekin/Etot[0], "b", label="Ekin")
ax1.plot(time, Eint/Etot[0], "g", label="Eint")
ax1.plot(time, Emag/Etot[0], "r", label="Emag")
ax1.plot(time, Etot/Etot[0], "k", label="Etot")
#ax2.plot(divB_err, "c", label="Error")

print(Etot)

ax1.set_xlabel("Time")
ax1.set_ylabel("Energy")
#ax2.set_ylabel("Divergence Error")

ax1.legend()

plt.savefig("EnergyPlot.png", dpi=200)
