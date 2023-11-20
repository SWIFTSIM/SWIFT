from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
import numpy as np
import sys
import os
import matplotlib
#matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.01)
    plt.colorbar(im, cax=cax)
    return

filename_base = "FastRotor_LR_"

nini=int(sys.argv[1])
nfin=int(sys.argv[2])
#for ii in range(61):
for ii in range(nini,nfin):

    print(ii)

    filename = filename_base + str(ii).zfill(4) + ".hdf5"
    data = load(filename)
    #print(data.metadata.gas_properties.field_names)
    boxsize = data.metadata.boxsize
    extent = [0, 1.0 , 0, 1.0]
    # cut the domian in half 
    data.metadata.boxsize*=[1.0,0.5,1.0]
    
    gas_gamma = data.metadata.gas_gamma
    print("Gas Gamma:",gas_gamma)
    if gas_gamma != 7./5.:
        print("WRONG GAS GAMMA")
        exit()

    mhdflavour = data.metadata.hydro_scheme["MHD Flavour"]
    mhd_scheme = data.metadata.hydro_scheme["MHD Scheme"]
    mhdeta = data.metadata.hydro_scheme["Resistive Eta"]
    git = data.metadata.code["Git Revision"]
    gitBranch = data.metadata.code["Git Branch"]
    scheme = data.metadata.hydro_scheme["Scheme"]
    kernel = data.metadata.hydro_scheme["Kernel function"]
    neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]
    
    try:
        dedhyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
        dedpar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
    except:
        dedhyp = 0.0
        dedpar = 0.0

    try:
        deddivV = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
        tensile = data.metadata.hydro_scheme["MHD Tensile Instability Correction Prefactor"]
        artdiff = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
    except:
        deddivV = 0.0
        artdiff = 0.0
        tensile = 1.0
 
    # First create a mass-weighted temperature dataset
    
    B = data.gas.magnetic_flux_densities
    divB = data.gas.magnetic_divergences
    P_mag = (B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2) / 2
    h = data.gas.smoothing_lengths
    
    normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
    
    data.gas.DivB_error = np.maximum(h * abs(divB) / normB, 1e-10)
    
    # Then create a mass-weighted B error dataset
    data.gas.mass_weighted_magnetic_divB_error = data.gas.masses * data.gas.DivB_error
    
    # Then create a mass-weighted B pressure dataset
    data.gas.mass_weighted_magnetic_pressures = data.gas.masses * P_mag

    # Then create a mass-weighted pressure dataset
    data.gas.mass_weighted_pressures = data.gas.masses * data.gas.pressures

    # Then create a mass-weighted Plasma Beta dataset
    data.gas.plasma_beta = data.gas.pressures / P_mag

    # Then create a mass-weighted speed dataset
    v = data.gas.velocities
    data.gas.mass_weighted_speeds = data.gas.masses * (
        v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2
    )
    # Then create a mass-weighted densities dataset
    data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities

    # Map in mass per area
    mass_map = slice_gas(data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024, project="masses", parallel=True  )
    
    # Map in density per area
    rho_map = slice_gas(
                 data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024, 
                 project="mass_weighted_densities", parallel=True)
    # Map in magnetism squared times mass per area
    mw_magnetic_pressure_map = slice_gas(
                 data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024,
                 project="mass_weighted_magnetic_pressures", parallel=True)

    # Map in pressure times mass per area
    mw_pressure_map = slice_gas(
                 data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024,
                 project="mass_weighted_pressures", parallel=True, )

    # Map in speed squared times mass per area
    mw_speeds_map = slice_gas(
                 data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024,
                 project="mass_weighted_speeds", parallel=True, )
    
    # Map in divB error times mass per area
    mw_ErrDivB_map = slice_gas(
                 data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024,
                 project="mass_weighted_magnetic_divB_error", parallel=True, )
    
    # Map in Plasma Beta times mass per area
    plasma_beta_map = slice_gas(
                 data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024,
                 project="plasma_beta", parallel=True, )

    rho_map = rho_map / mass_map
    magnetic_pressure_map = mw_magnetic_pressure_map / mass_map
    speed_map = mw_speeds_map / mass_map
    pressure_map = mw_pressure_map / mass_map
    ErrDivB_map = mw_ErrDivB_map / mass_map
    #plasma_beta_map

    #fig = plt.figure(figsize=(12, 11), dpi=100)
    fig = plt.figure(figsize=(12, 8), dpi=100)
    
    ax1 = fig.add_subplot(231)
    im1 = ax1.imshow(rho_map.T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm(vmax=14,vmin=1))
    ax1.set_title("Density")
    set_colorbar(ax1, im1)
    
    ax2 = fig.add_subplot(232)
    im2 = ax2.imshow(magnetic_pressure_map.T, origin="lower", extent=extent, cmap="plasma", norm=Normalize(vmax=3,vmin=0.1))
    ax2.set_title("Magnetic Pressure")
    set_colorbar(ax2, im2)
    
    ax3 = fig.add_subplot(233)
    im3 = ax3.imshow(speed_map.T, origin="lower", extent=extent, cmap="cividis", norm=Normalize(vmax=1.6,vmin=0.1))
    ax3.set_title("v^2")
    set_colorbar(ax3, im3)
    
    ax4 = fig.add_subplot(234)
    im4 = ax4.imshow(pressure_map.T, origin="lower", extent=extent, cmap="viridis", norm=Normalize(vmax=2.1,vmin=0.1))
    ax4.set_title("Internal Pressure")
    set_colorbar(ax4, im4)
    
    #ax5 = fig.add_subplot(235)
    #im5 = ax5.imshow(plasma_beta_map.T, origin="lower", extent=extent, cmap="magma", norm=LogNorm())
    #ax5.set_title("plasma beta")
    #set_colorbar(ax5, im5)
    
    ax5 = fig.add_subplot(235)
    im5 = ax5.imshow(ErrDivB_map.T, origin="lower", extent=extent, cmap="gray", norm=LogNorm(vmax=1,vmin=0.001))
    ax5.set_title("Err(DivB)")
    set_colorbar(ax5, im5)

    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xlabel("x ")
        ax.set_ylabel("y ")
    
    ax6 = fig.add_subplot(236)
    
    text_fontsize = 8
    ax6.text(
        0.1,
        0.9,
        "Fast Rotor $t=%.2f$" % data.metadata.time,
        fontsize=text_fontsize,
    )
    ax6.text(0.1, 0.85, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
    ax6.text(0.1, 0.8, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize)
    ax6.text(0.1, 0.75, scheme.decode("utf-8"), fontsize=text_fontsize)
    ax6.text(0.1, 0.7, kernel.decode("utf-8"), fontsize=text_fontsize)
    ax6.text(0.1, 0.65, "$%.2f$ neighbours" % (neighbours), fontsize=text_fontsize)
    ax6.text(
        0.1,
        0.55,
        "$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:25],
        fontsize=text_fontsize,
    )
    ax6.text(0.1, 0.5, "$Resitivity_\\eta:%.4f$ " % (mhdeta), fontsize=text_fontsize)
    ax6.text(0.1, 0.45, "$Dedner Parameters: $", fontsize=text_fontsize)
    ax6.text(0.1, 0.4, "$[hyp, par, div] [:%.3f,%.3f,%.3f]$ " % (dedhyp,dedpar,deddivV), fontsize=text_fontsize)
    ax6.text(0.1, 0.35, "$Tensile Prefactor:%.4f$ " % (tensile), fontsize=text_fontsize)
    ax6.text(0.1, 0.3, "$Art. Diffusion:%.4f$ " % (artdiff), fontsize=text_fontsize)
    ax6.tick_params(left=False, right=False, labelleft=False, labelbottom=False, bottom=False)
    #ax6.set_xlabel("")
    #ax6.plot(frameon=False)

    plt.tight_layout()
    plt.savefig(filename_base + str(ii).zfill(4) + ".jpg")
    plt.close()
