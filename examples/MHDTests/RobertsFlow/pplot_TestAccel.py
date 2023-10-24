from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas

import numpy as np
import matplotlib.pyplot as plt
import sys

def plotx4(f_in,f_out):
   data = load(f_in)
   print(f_in)
   #print(f_in+" "+data.metadata.time)
   #print(data.metadata)
   #print(data.metadata.gas_properties.field_names)
   #quit()

   mhdflavour = data.metadata.hydro_scheme["MHD Flavour"]
   #dedhyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
   #dedpar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
   #mhdeta = data.metadata.hydro_scheme["Resistive Eta"]
   git = data.metadata.code["Git Revision"]
   gitBranch = data.metadata.code["Git Branch"]
   scheme = data.metadata.hydro_scheme["Scheme"]
   kernel = data.metadata.hydro_scheme["Kernel function"]
   neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]

   r_gas = data.gas.coordinates
   rho_gas = data.gas.densities
   B_gas = data.gas.magnetic_flux_density
   h_gas = data.gas.smoothing_lengths
   divB_gas = data.gas.magnetic_divergence
   vel_gas = data.gas.velocities

   x_gas = r_gas[:, 0]
   y_gas = r_gas[:, 1]

   B_mag_gas = (
      B_gas[:, 0] * B_gas[:, 0] + B_gas[:, 1] * B_gas[:, 1] + B_gas[:, 2] * B_gas[:, 2]
   )
   B_mag_gas = B_mag_gas / 2

   vel_mag_gas = (
       vel_gas[:, 0] * vel_gas[:, 0] + vel_gas[:, 1] * vel_gas[:, 1] + vel_gas[:, 2] * vel_gas[:, 2]
   )

# N = 100000

   x_val = x_gas.value  # [:N]
   y_val = y_gas.value  # [:N]

   rho_val = rho_gas.value  # [:N]

   B_mag_val = np.where(B_mag_gas.value != 0.0, B_mag_gas.value, 1E-9)

   V_mag_val = vel_mag_gas.value

   h_val = h_gas.value
   divB_val = np.where(divB_gas.value != 0.0, divB_gas.value, 1E-18)


   err_val = np.log10(h_val * abs(divB_val) / np.sqrt(2 * B_mag_val) + 1e-8)

# err_val[np.isnan(err_val)] = 0

   fig, ax = plt.subplots(1, 4, figsize=(16, 5))
#######################
   cs1 = ax[0].tricontourf(
      x_val, y_val, rho_val, levels=np.linspace(0.0, 5.0, 200), cmap="gist_heat"
   )
   cbar1 = plt.colorbar(
      cs1, ax=ax[0], ticks=np.linspace(0.0, 5.0, 10), fraction=0.046, pad=0.04
   )
   cbar1.set_label(r"$\rho$")
   ax[0].set(xticks=[], yticks=[], aspect="equal")
#####################
   cs2 = ax[1].tricontourf(
      #x_val, y_val, B_mag_val, levels=np.linspace(0.0, 0.4, 200), cmap="turbo"
      x_val, y_val, B_mag_val, cmap="plasma"
   )
   cbar2 = plt.colorbar(
    #cs2, ax=ax[1], ticks=np.linspace(0.0, 0.4, 5), fraction=0.046, pad=0.04
    cs2, ax=ax[1], fraction=0.046, pad=0.04
   )
   cbar2.set_label(r"$B^2 / 2$")
############
   ax[1].set(xticks=[], yticks=[], aspect="equal")

   cs3 = ax[2].tricontourf(
    x_val, y_val, err_val, levels=np.linspace(-4.0, 1.0, 200), cmap="gray"
   )
   cbar3 = plt.colorbar(
    cs3, ax=ax[2], ticks=np.linspace(-4.0, 1.0, 6), fraction=0.046, pad=0.04
   )
   cbar3.set_label(r"$Error$")
   ax[2].set(xticks=[], yticks=[], aspect="equal")
#######
   cs4 = ax[3].tricontourf(
    #x_val, y_val, V_mag_val, levels=np.linspace(0,1.5,200),cmap="viridis"
    x_val, y_val, V_mag_val,cmap="viridis"
   )
   cbar4 = plt.colorbar(
    #cs4, ax=ax[3], ticks=np.linspace(0.0, 1.5, 7), fraction=0.046, pad=0.04
    cs4, ax=ax[3], fraction=0.046, pad=0.04
   )
   cbar4.set_label(r"$Velocity$")
   ax[3].set(xticks=[], yticks=[], aspect="equal")
########

   fig.subplots_adjust(hspace = 0.02)
   fig.subplots_adjust(wspace = 0.4)
   fig.tight_layout()

   text_fontsize = 7
   ax[0].text(0.1, 1.25, "Orzag Tang Vortex Time $t=%.2f$" % data.metadata.time, fontsize=text_fontsize)
   ax[0].text(0.1, 1.2, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
   ax[0].text(0.1, 1.15, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize)
   ax[0].text(0.1, 1.1, scheme.decode("utf-8"), fontsize=text_fontsize)
   ax[0].text(0.1, 1.05, kernel.decode("utf-8"), fontsize=text_fontsize)
   ax[0].text(0.6, 1.05, "$%.2f$ neighbours" % (neighbours), fontsize=text_fontsize)
   ax[1].text(0.1,1.15,"$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:30],fontsize=text_fontsize)
   #ax[1].text(0.1, 1.1, "$Resitivity_\\eta:%.4f$ " % (mhdeta), fontsize=text_fontsize)

   plt.savefig(f_out)
   plt.close()


n_ini=000
n_fin=50


fname_in  = [sys.argv[1]+'RobertsFlow_{:04}.hdf5'.format(i) for i in range(n_ini, n_fin)]
fname_out = [sys.argv[1]+'RF_{:04}.jpg'.format(i) for i in range(n_ini, n_fin)]

for i in range(n_ini,n_fin):
   plotx4(fname_in[i],fname_out[i])

