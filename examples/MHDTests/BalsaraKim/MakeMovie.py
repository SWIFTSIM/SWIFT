from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
from swiftsimio.visualisation.slice import slice_gas

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, BoundaryNorm, Normalize
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch

import unyt
import glob
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser(description='make a move of the snapshots')
parser.add_argument('snapdir', type=str,
                    help='the directory containing the snapshots')
parser.add_argument('z', type=float,
                    help='the relative slice height')
parser.add_argument('savefile', type=str,
                    help='the name of the savefile')
parser.add_argument('--snapbasename', type=str, nargs='?', default='BK',
                    help='the basename of the snapshots')
args = parser.parse_args()

img_res = 128

def abs_vec(vec):
    if len(vec.shape)==2:
        res = np.sqrt(vec[:, 0] ** 2 + vec[:, 1] ** 2 + vec[:, 2] ** 2)
    elif len(vec.shape)==3:
        res = np.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2)
    return res

def rms_vec(vec):
    res = np.sqrt(
        np.mean(vec[:, 0] * vec[:, 0] + vec[:, 1] * vec[:, 1] + vec[:, 2] * vec[:, 2])
    )
    return res

def set_colorbar(fig, ax, im):
    """sets a nicely placed colorbar at the right of the given axis object <ax>, in figure <fig>
        with image <im>"""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, extend='both')
    cbar.ax.tick_params(labelsize=12)
    return

def convert_to_gauss (B, micro=True):
    """convert magnetic flux density cosmo array to units gauss"""
    B.convert_to_units(unyt.tesla)
    B.convert_to_units(unyt.microgauss)
    return B

def convert_to_gauss_pc (divB, micro=True):
    """convert magnetic divergence cosmo array to units gauss.
        unyt can't handle this conversion."""
    divB.convert_to_units(unyt.T / unyt.pc)
    divB *= 1*unyt.pc

    divB.convert_to_units(unyt.microgauss)

    divB /= 1*unyt.pc
    return divB

def make_slice(key, z):
    res = slice_gas(data,
                    z_slice=z,
                    resolution=img_res,
                    project=key,
                    parallel=True,
                    periodic=True)
    return res

def prepare_sliced_quantity(quantity, z, mass_map, isvec=False):
    mass = data.gas.masses
    if isvec:
        # Prepare vector quantity for slicing
        data.gas.mass_weighted_temp_qx = mass * quantity[:, 0]
        data.gas.mass_weighted_temp_qy = mass * quantity[:, 1]
        data.gas.mass_weighted_temp_qz = mass * quantity[:, 2]

        # Apply slicing for each component
        sliced_quantity_x = make_slice("mass_weighted_temp_qx", z=z) / mass_map
        sliced_quantity_y = make_slice("mass_weighted_temp_qy", z=z) / mass_map
        sliced_quantity_z = make_slice("mass_weighted_temp_qz", z=z) / mass_map

        # Join components together (should be optimized!!!)
        the_units = sliced_quantity_x.units
        sliced_quantity = np.array([sliced_quantity_x.value,
                                    sliced_quantity_y.value,
                                    sliced_quantity_z.value])
        sliced_quantity = sliced_quantity * the_units
    else:
        # prepare scalar quantity for slicing
        data.gas.mass_weighted_temp_q = data.gas.masses * quantity[:]
        # apply slicing
        sliced_quantity = make_slice("mass_weighted_temp_q", z=z) / mass_map
        # Convert scalar quantity to physical
        sliced_quantity = sliced_quantity
    return sliced_quantity

def make_projection (key):
    res = project_gas(data,
                      resolution=img_res,
                      project=key,
                      parallel=True,
                      periodic=True)
    return res

def prepare_projected_quantity(quantity, mass_map, isvec=False):
    mass = data.gas.masses
    if isvec:
        # Prepare vector quantity for slicing
        data.gas.mass_weighted_temp_qx = mass * quantity[:, 0]
        data.gas.mass_weighted_temp_qy = mass * quantity[:, 1]
        data.gas.mass_weighted_temp_qz = mass * quantity[:, 2]

        # Apply slicing for each component
        sliced_quantity_x = make_projection("mass_weighted_temp_qx") / mass_map
        sliced_quantity_y = make_projection("mass_weighted_temp_qy") / mass_map
        sliced_quantity_z = make_projection("mass_weighted_temp_qz") / mass_map

        # Convert vector quantity to physical
        sliced_quantity_x = sliced_quantity_x.to_physical()
        sliced_quantity_y = sliced_quantity_y.to_physical()
        sliced_quantity_z = sliced_quantity_z.to_physical()

        # Join components together (should be optimized!!!)
        the_units = sliced_quantity_x.units
        sliced_quantity = np.array([sliced_quantity_x.value,
                                    sliced_quantity_y.value,
                                    sliced_quantity_z.value])
        sliced_quantity = sliced_quantity * the_units
    else:
        # prepare scalar quantity for slicing
        data.gas.mass_weighted_temp_q = data.gas.masses * quantity[:]
        # apply slicing
        sliced_quantity = make_projection("mass_weighted_temp_q") / mass_map
        # Convert scalar quantity to physical
        sliced_quantity = sliced_quantity.to_physical()
    return sliced_quantity

### INPUT ###
snapdir = args.snapdir
z_relative = args.z

#initialise lists
rho_projs  = []

rho_slices = []

Babs_slices = []
B_slices    = []

vabs_slices = []
v_slices    = []

R0_slices = []

t = []

#get filenames
fnames = sorted(glob.glob(snapdir + args.snapbasename + '_????.hdf5'))

#looping over all snapshots to get the data for the movie
for i in tqdm(range(len(fnames)), desc='computing slices'):
    data = load(fnames[i])

    boxsize = data.metadata.boxsize

    #obtaining the quantities
    time = data.metadata.time
    time.convert_to_units(unyt.Myr)
    t.append(np.round(time.value, 3))

    rho = data.gas.densities.to_physical() / unyt.mh
    rho.convert_to_units(unyt.cm**-3)

    B = convert_to_gauss(data.gas.magnetic_flux_densities.to_physical())

    v = data.gas.velocities.to_physical()

    divB = convert_to_gauss_pc(data.gas.magnetic_divergences.to_physical())
    h    = data.gas.smoothing_lengths.to_physical()
    R0   = np.abs( h*divB / abs_vec(B) ).value

    #slicing preperation
    masses = make_slice('masses', z=z_relative*boxsize[2])

    #slicing quantities
    rho_sl = prepare_sliced_quantity(rho, z=z_relative*boxsize[2], mass_map=masses)
    R0 = prepare_sliced_quantity(R0, z=z_relative*boxsize[2], mass_map=masses)

    B = prepare_sliced_quantity(B, z=z_relative*boxsize[2], mass_map=masses, isvec=True)
    v = prepare_sliced_quantity(v, z=z_relative*boxsize[2], mass_map=masses, isvec=True)

    #projection preperation
    masses = make_projection("masses")
    masses.convert_to_units(1e8 * unyt.Msun / unyt.kpc**2)

    #storing the slices
    rho_slices.append(rho_sl)

    Babs_slices.append(abs_vec(B))
    B_slices.append(B)

    v.convert_to_units(unyt.km/unyt.s)
    vabs_slices.append(abs_vec(v))
    v_slices.append(v)

    R0_slices.append(R0)

    #storing the projections
    rho_projs.append(masses)

# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"][0]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"][0]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"][0]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"][0]
eta = data.metadata.hydro_scheme["Resistive Eta"][0]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"][0]
Np = data.metadata.n_gas

#good colormap scales (finding them (min,max) overloads kernel)
rho_vmin, rho_vmax = 5e-1, 5e2
B_vmin, B_vmax = 1, 1e2
v_vmin, v_vmax = 1, 1e3

rho_pr_vmin, rho_pr_vmax = 1.7e-1, 2.7

#bounds for the R0 plot colorbar
bounds = np.linspace(-3, 1, 9)
cbar_ticks = np.linspace(-3, 1, 5)

#defining the grid
dimy = len(masses)
dimx = len(masses[0])
new_x = np.linspace(0, boxsize[0].value, dimx)
new_y = np.linspace(0, boxsize[1].value, dimy)

extent = np.array([0, boxsize[0].value,
                   0, boxsize[1].value])

### starting the animation ###

plt.rcParams.update({
    "axes.labelsize": 14,
    "axes.titlesize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10
})

fig, axes = plt.subplots(3, 2, figsize=(10, 15), sharex=True, sharey=True)
fig.subplots_adjust(wspace=.25)

text_common_args = dict(fontsize=10, ha="center", va="center", transform=axes[0,1].transAxes)
time_text = axes[0, 1].text(0.5, 0.8, '', **text_common_args)   # outside init to access time_text in update function

bounds = np.linspace(-3, 1, 9)
cbar_ticks = np.linspace(-3, 1, 5)

def init ():
    #setting the axis limits and labels
    for ax in axes.flatten():
        ax.set_xlim(0, extent[1])
        ax.set_ylim(0, extent[3])

    plt.setp(axes[-1,:], xlabel=r'$x \ \mathrm{[kpc]}$')
    plt.setp(axes[:,0], ylabel=r'$y \ \mathrm{[kpc]}$')

    # initialise text
    axes[0, 1].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
    axes[0, 1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
    axes[0, 1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
    axes[0, 1].text(0.5,0.4,kernel.decode("utf-8") + " with %.2f neighbours" % (neighbours),
        **text_common_args)
    axes[0, 1].text(0.5, 0.3, "Artificial diffusion: %.2f " % (artDiffusion), **text_common_args)
    axes[0, 1].text(0.5,0.2,"Dedner Hyp, Hyp_div(v), Par: %.2f,%.2f,%.2f " % (dedHyp, dedHypDivv, dedPar),
        **text_common_args)
    axes[0, 1].text(0.5, 0.1, r'Physical resistivity $\eta$: %.2f ' % (eta), **text_common_args)
    axes[0, 1].text(0.5, 0.0, r'Number of particles $N_p$: %d ' % (Np), **text_common_args)
    axes[0, 1].axis("off")

    #placing the colorbars, will be the same for the whole movie
    set_colorbar(fig, axes[0,0], ScalarMappable(norm=Normalize(vmin=rho_pr_vmin, vmax=rho_pr_vmax), cmap='viridis'))
    set_colorbar(fig, axes[1,0], ScalarMappable(norm=LogNorm(rho_vmin, rho_vmax), cmap='viridis'))
    set_colorbar(fig, axes[1,1], ScalarMappable(norm=LogNorm(v_vmin, v_vmax), cmap='inferno'))
    set_colorbar(fig, axes[2,0], ScalarMappable(norm=LogNorm(B_vmin, B_vmax), cmap='magma'))
    set_colorbar(fig, axes[2,1], ScalarMappable(norm=BoundaryNorm(boundaries=bounds, ncolors=256,
                                               extend='both'), cmap='jet'))

    #setting the titles
    axes[0,0].set_title(r'$\Sigma \ \mathrm{[10^8 \ M_\odot \ kpc^{-2}]}$')
    axes[1,0].set_title(r'$\rho/m_H \ \mathrm{[cm^{-3}]}$')
    axes[1,1].set_title(r'$v \ \mathrm{[km / s]}$')
    axes[2,0].set_title(r'$B \ \mathrm{[\mu G]}$')
    axes[2,1].set_title(r'$h|\nabla \cdot \mathbf{B}|/B$')

    return

def update (frame):
    #removing the previous streamplots, otherwise streamplots will be drawn on top of each other
    if frame > 0:   #condition to make sure there is a streamplot to remove
        for ax in [axes[1,1], axes[2,0]]:
            ax.collections[0].remove() #removing the lines

            for artist in ax.get_children():
                if isinstance(artist, FancyArrowPatch):
                    artist.remove() #removing the arrowheads

    # update time 
    time_text.set_text(r'Balsara-Kim test at time $t=%.2f$ Myr' % t[frame])

    #drawing the images
    axes[0,0].imshow(rho_projs[frame].value.T, vmin=rho_pr_vmin, vmax=rho_pr_vmax,
                     extent=extent, origin='lower', cmap='viridis')
    axes[1,0].imshow(rho_slices[frame].value.T, norm=LogNorm(rho_vmin, rho_vmax),
                     extent=extent, origin='lower', cmap='viridis')
    axes[1,1].imshow(vabs_slices[frame].value.T, norm=LogNorm(vmin=v_vmin, vmax=v_vmax),
                     extent=extent, origin='lower', cmap='inferno')
    axes[2,0].imshow(Babs_slices[frame].value.T, norm=LogNorm(B_vmin, B_vmax),
                     extent=extent, origin='lower', cmap='magma')
    axes[2,1].imshow(np.log10(R0_slices[frame].value.T), 
                     norm=BoundaryNorm(boundaries=bounds, ncolors=256, extend='both'),
                     extent=extent, origin='lower', cmap='jet')

    #drawing the field lines
    axes[1,1].streamplot(new_x, new_y,
                         v_slices[frame][0].value.T, v_slices[frame][1].value.T,
                         color='w', density=1.5, linewidth=.5, arrowsize=.5)
    axes[2,0].streamplot(new_x, new_y,
                         B_slices[frame][0].value.T, B_slices[frame][1].value.T,
                         color='w', density=1.5, linewidth=.5, arrowsize=.5)

    return

#executing the movie
ani = FuncAnimation(fig, update, frames=tqdm(range(len(fnames)), desc='animating movie'),
                    init_func=init, blit=False, interval=250)

ani.save(args.savefile)





