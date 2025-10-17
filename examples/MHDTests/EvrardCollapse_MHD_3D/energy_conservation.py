from swiftsimio import load_statistics
import matplotlib.pyplot as plt
import sys

data = load_statistics("statistics_vsig_dt.txt")

# print(data)

N = 300

Ekin = data.kin_energy.value[:N]
Eint = data.int_energy.value[:N]
Emag = data.mag_energy.value[:N]
Epot = data.pot_energy.value[:N]

Etot = Ekin + Eint + Emag + Epot

Enorm = 1.0  # Etot[0]

divB_err = data.divb_err.value[:N]

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

# ax[0].plot(Ekin/Enorm, 'b', label='Ekin')
# ax[0].plot(Eint/Enorm, 'g', label='Eint')
ax[0].plot(Emag / Enorm, "r", label="Emag")
# ax[0].plot(Epot/Enorm, 'm', label='Epot')
# ax[0].plot(Etot/Enorm, 'k', label='Etot')
ax[1].plot(divB_err / 10000, "c", label="Error")

ax[0].set_xlabel("Time")
ax[0].set_ylabel("Energy")
ax[1].set_xlabel("Time")
ax[1].set_ylabel("Divergence Error")

ax[0].grid(alpha=0.2)
ax[1].grid(alpha=0.2)

ax[0].legend()
ax[1].legend()

plt.tight_layout()

plt.savefig(sys.argv[1], dpi=200)
