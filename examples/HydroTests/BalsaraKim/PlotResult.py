import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas

import unyt

plt.rcParams.update({
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 20,
    "savefig.dpi": 100
})

img_res = 512
div_reg = 1e-30

def set_colorbar(fig, ax, im):
    """sets a nicely placed colorbar at the right of the given axis object <ax>, in figure <fig>
        with image <im>"""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, extend='both')
    cbar.ax.tick_params(labelsize=12)
    return

def make_projection (data, key):
    res = project_gas(data,
                      resolution=img_res,
                      project=key,
                      parallel=True,
                      periodic=True)
    return res

def make_column_density (fname):      
    data = load(fname)

    # get some metadata 
    global boxsize
    boxsize = data.metadata.boxsize
    time = data.metadata.time
    time.convert_to_units(unyt.Myr)

    # make projection 
    masses = make_projection(data, "masses")
    
    extent = np.array([0, boxsize[0].value,
                       0, boxsize[1].value])

    return time, masses, extent

# the snapshots to show
snapnrs = [10, 20, 38]
fnames = ['snapshots/BK_' + str(i).zfill(4) + '.hdf5' for i in snapnrs]

# prepare and store the column densities
times = []
col_densities = []

for fname in fnames:
    t, rho, extent = make_column_density(fname)
    times.append(t)
    col_densities.append(rho)

# parameters for scalebar to plot
L_box = boxsize[0]
L_box.convert_to_units(unyt.pc)

perc = 0.25     # percentage of the boxsize to use as scale

### PLOTTING ###

fig, axes = plt.subplots(1, 3, figsize=(12,4))

for i in range(3):
    im_rho = axes[i].imshow(col_densities[i].value.T/1e8, vmin=1.7e-1, vmax=2.7, 
                            extent=extent, origin='lower', cmap='viridis')

    axes[i].set_title(r'$t=$'+str(np.round(times[i],3)) + ' Myr')

    # make nice axes edges
    for loc in ['bottom', 'top', 'right', 'left']:
        axes[i].spines[loc].set_color('white')
    axes[i].tick_params(which='both', axis='both', direction='in', color='white', top=True, right=True)

    axes[i].set_yticks(np.linspace(0, 0.2, 9))
    axes[i].set_xticks(np.linspace(0, 0.2, 9))
    
    axes[i].set_xticklabels([])
    axes[i].set_yticklabels([])
    
    # make scalebar
    axes[i].hlines(0.05, 0.75-perc/2, 0.75+perc/2,
                   transform=axes[i].transAxes, color='white')
    axes[i].text(0.75 , 0.06, str(int(perc*L_box.value)) + ' pc',
                 transform=axes[i].transAxes, 
                 color='white', horizontalalignment='center',
                 verticalalignment='bottom', size=11)

# make colorbar
cbar_ax = fig.add_axes([0.905, 0.028, 0.02, 0.895])
cbar = fig.colorbar(im_rho, cax=cbar_ax,
                    label=r'$\Sigma \ \mathrm{[10^8 \ M_\odot  kpc^{-2}]}$',
                    format=ScalarFormatter())
cbar.minorformatter = ScalarFormatter()    

# adjust spacing between axes
fig.subplots_adjust(left=0, right=0.9, top=0.95, bottom=0., wspace=0, hspace=0)

# SAVING...
fig.savefig('BalsaraKim-ColumnDensity.jpg')
