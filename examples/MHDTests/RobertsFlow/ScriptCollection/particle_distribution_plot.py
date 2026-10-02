import numpy as np
#import h5py
import unyt
import matplotlib.pyplot as plt
import argparse
import scipy as sp

from swiftsimio import load
from swiftsimio.visualisation.volume_render import render_gas
from swiftsimio.visualisation.slice import slice_gas
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",    # Computer Modern
})

if __name__ == "__main__":

    import argparse as ap

    parser = ap.ArgumentParser(description="Generate ICs for ABC flow")

    parser.add_argument(
        "-in",
        "--input_filename",
        help="source file name",
        default="",
        type=str,
    )
    parser.add_argument(
        "-outim",
        "--out_image",
        help="output image name",
        default="im.png",
        type=str,
    )
    parser.add_argument(
        "-z",
        "--z_slab",
        help="slab center position",
        default=0.5,
        type=float,
    )
    parser.add_argument(
        "-d",
        "--d_slab",
        help="slab height",
        default=None,
        type=float,
    )

    args = parser.parse_args()
 
    data = load(args.input_filename)

    coordinates = data.gas.coordinates
    h = data.gas.smoothing_lengths
    rho = data.gas.densities
    m = data.gas.masses
   
    h_mean = np.mean(h)
    d_ip_mean = np.mean((m/rho)**(1/3))
    Lbox = data.metadata.boxsize

    z_cut = args.z_slab * coordinates.units
    d_slab = args.d_slab

    if d_slab==None:
        d_cut = d_ip_mean
    else:
        d_cut = args.d_slab * coordinates.units

    if ( (z_cut > d_cut/2.0) | (Lbox[2]-z_cut<d_cut/2.0)):
        mask = ( (coordinates[:,2]>=z_cut-d_cut/2.0) & (coordinates[:,2]<=z_cut+d_cut/2.0) )
    elif (z_cut < d_cut/2.0):
        mask_bottom = ( (coordinates[:,2]>=0.0) & (coordinates[:,2]<=z_cut+d_cut/2.0) ) 
        mask_top = ( (coordinates[:,2]>=Lbox[2] - (d_cut/2.0-z_cut) ) & (coordinates[:,2]<=Lbox[2]) ) 
        mask = ( mask_bottom & mask_top )
    elif (Lbox-z_cut > d_cut/2.0):
        mask_bottom = ( (coordinates[:,2]>=z_cut - d_cut/2.0 ) & (coordinates[:,2]<=Lbox[2]) ) 
        mask_top = ( (coordinates[:,2]>=0.0 ) & (coordinates[:,2]<=d_cut/2.0 + (-Lbox[2]+z_cut)) ) 
        mask = ( mask_bottom & mask_top )
    else:
        print('Error in slicing')

    coordinates = coordinates[mask]
    rho = rho[mask]
    drho_over_rho = 1 - np.mean(rho)/rho

    # Generate figure
    nx = 1
    ny = 1

    lx = 10#6
    ly = 8.0#5

    fig, ax = plt.subplots(ny, nx, sharey=True,sharex=True, figsize=(lx * nx, ly * ny))

    from matplotlib.ticker import AutoMinorLocator
    # Ticks on all four sides
    ax.tick_params(axis='both',
                  which='major',
                  direction='in',
                  top=True,
                  right=True,
                  length=7,
                  width=1.2,
                  labelsize=16)
    
    ax.tick_params(axis='both',
                  which='minor',
                  direction='in',
                  top=True,
                  right=True,
                  length=4,
                  width=1.0)
    
    # Add minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
 
    x = coordinates[:,0]
    y = coordinates[:,1]

    sc = ax.scatter(
        x, y,
        #s=0.25,
        s=3.0, #0.5,
        c=drho_over_rho,
        cmap='bwr',
        vmin=-0.15,
        vmax=0.15
    )
    cbar = fig.colorbar(sc, ax=ax)
    #cbar.set_label(r'$\delta \rho / \rho$', fontsize=15)
    cbar.set_ticks([-0.15,-0.1,-0.05, 0.0,0.05,0.1,0.15])
    cbar.ax.tick_params(
        labelsize=15,  # font size of tick labels
        #length=6,      # tick length
        #width=1.5      # tick width
    )

    # Add ticks and labels
    locs_y = [Lbox[0]/4 , Lbox[0]/2, 3 * Lbox[0]/4, Lbox[0]]
    locs_x = [0.0 ,Lbox[0]/4 , Lbox[0]/2, 3 * Lbox[0]/4, Lbox[0]]
    labels_y = [r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$',r'$2\pi$']
    labels_x = [r'$ 0 $',r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$',r'$2\pi$']
    ax.set_yticks(locs_y, labels_y,fontsize=25)
    ax.set_xticks(locs_x, labels_x,fontsize=25)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(0, Lbox[0])
    ax.set_ylim(0, Lbox[1])
    #ax.set_title(r"$\delta \rho/\rho \;$",fontsize=30,pad=10)

    fig.tight_layout()
    plt.savefig(args.out_image, dpi=200)

