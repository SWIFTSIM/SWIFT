from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas

import numpy as np
import matplotlib.pyplot as plt
import sys

data = load(sys.argv[1])
print(data.metadata.time)
# print(data.metadata)
# quit()
# print(data.metadata.gas_properties.field_names)

mhdflavour = data.metadata.hydro_scheme["MHD Flavour"]
# dedhyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
# dedpar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
mhdeta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
scheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]


r_gas = data.gas.coordinates
rho_gas = data.gas.densities
B_gas = data.gas.magnetic_flux_densities
h_gas = data.gas.smoothing_lengths
divB_gas = data.gas.magnetic_divergences

x_gas = r_gas[:, 0]
y_gas = r_gas[:, 1]

B_mag_gas = (
    B_gas[:, 0] * B_gas[:, 0] + B_gas[:, 1] * B_gas[:, 1] + B_gas[:, 2] * B_gas[:, 2]
)
B_mag_gas = B_mag_gas / 2

# N = 100000

x_val = x_gas.value  # [:N]
y_val = y_gas.value  # [:N]

rho_val = rho_gas.value  # [:N]

B_mag_val = B_mag_gas.value  # [:N]

h_val = h_gas.value
divB_val = divB_gas.value

err_val = np.log10(h_val * abs(divB_val) / np.sqrt(2 * B_mag_val) + 1e-8)

# err_val[np.isnan(err_val)] = 0

fig, ax = plt.subplots(1, 3, figsize=(16, 5))

cs1 = ax[0].tricontourf(
    x_val, y_val, rho_val, levels=np.linspace(0.0, 0.5, 200), cmap="gist_heat"
)
cbar1 = plt.colorbar(
    cs1, ax=ax[0], ticks=np.linspace(0.0, 0.5, 6), fraction=0.046, pad=0.04
)
cbar1.set_label(r"$\rho$")
ax[0].set(xticks=[], yticks=[], aspect="equal")

cs2 = ax[1].tricontourf(
    x_val, y_val, B_mag_val, levels=np.linspace(0.0, 0.4, 200), cmap="turbo"
)
cbar2 = plt.colorbar(
    cs2, ax=ax[1], ticks=np.linspace(0.0, 0.4, 5), fraction=0.046, pad=0.04
)
cbar2.set_label(r"$B^2 / 2$")
ax[1].set(xticks=[], yticks=[], aspect="equal")

cs3 = ax[2].tricontourf(
    x_val, y_val, err_val, levels=np.linspace(-4.0, 1.0, 200), cmap="gray"
)
cbar3 = plt.colorbar(
    cs3, ax=ax[2], ticks=np.linspace(-4.0, 1.0, 6), fraction=0.046, pad=0.04
)
cbar3.set_label(r"$Error$")
ax[2].set(xticks=[], yticks=[], aspect="equal")

# fig.subplots_adjust(hspace = 0.0)
fig.subplots_adjust(wspace=0.4)

text_fontsize = 7
ax[0].text(
    0.1,
    1.25,
    "Orzag Tang Vortex Time $t=%.2f$" % data.metadata.time,
    fontsize=text_fontsize,
)
ax[0].text(0.1, 1.2, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
ax[0].text(0.1, 1.15, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize)
ax[0].text(0.1, 1.1, scheme.decode("utf-8"), fontsize=text_fontsize)
ax[0].text(0.1, 1.05, kernel.decode("utf-8"), fontsize=text_fontsize)
ax[0].text(0.6, 1.05, "$%.2f$ neighbours" % (neighbours), fontsize=text_fontsize)
ax[1].text(
    0.1,
    1.15,
    "$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:30],
    fontsize=text_fontsize,
)
ax[1].text(0.1, 1.1, "$Resitivity_\\eta:%.4f$ " % (mhdeta), fontsize=text_fontsize)

plt.savefig(sys.argv[2])
quit()

# monopole_map = project_gas(data, resolution=512, project="monopole_term", parallel=True)
B_map = project_gas(
    data, resolution=512, project="magnetic_flux_densities", parallel=True
)


from matplotlib.pyplot import imsave

# from matplotlib.colors import LogNorm

imsave("test2.png", B_map.value, cmap="Spectral", vmin=0.0, vmax=0.15)

quit()
