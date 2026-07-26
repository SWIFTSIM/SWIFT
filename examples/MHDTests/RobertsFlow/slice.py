import numpy as np
#import h5py
import unyt
import matplotlib.pyplot as plt
import argparse
import scipy as sp

from swiftsimio import load
#from swiftsimio.visualisation.volume_render import render_gas
from swiftsimio.visualisation.slice import slice_gas
from matplotlib.image import imread
from matplotlib.colors import LinearSegmentedColormap


# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
#argparser.add_argument("resolution")
args = argparser.parse_args()



# Load snapshot
filename = args.input
data = load(filename)

time = np.round(data.metadata.time,2)
print("Plotting at ", time)

# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
eta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]

Lbox = data.metadata.boxsize

# creating grid for B field
res = 64

def periodic_continue(Q):

    """
    Periodically continue the image
    Parameters: Q: quantity, with size (N_x, N_y)
    Returns: 
                Qp: quantity, with size (N_x+1, N_y+1)
    """

    M, N = Q.shape
    #Q = np.arange(M*N).reshape(M, N)  # example data

    # Periodic continuation
    Qp = np.zeros((M+1, N+1), dtype=Q.dtype)
    Qp[:-1, :-1] = Q         # copy original
    Qp[-1, :-1] = Q[0, :]     # wrap last row
    Qp[:-1, -1] = Q[:, 0]     # wrap last column
    Qp[-1, -1]  = Q[0, 0]     # wrap bottom-right corner
    return Qp

def make_2D_slice_z(data,
                    height
                    ):

    """
    Making SPH B slice of the box at specific height
    Parameters: data: SPH snapshot data
                height: slice height
    Returns: 
                dict: B - magnetic field array of size (N_x, N_y, 3), z - slice height
    """

    visualise_region_xy = [
    0*unyt.cm,
    Lbox[0],
    0*unyt.cm,
    Lbox[1],
    ]

    common_arguments_xy = dict(
        data=data,
        resolution=res,
        parallel=True,
        region=visualise_region_xy,
        z_slice=height,  
        periodic=True,
    )

    mass_map_ij  = slice_gas(**common_arguments_xy, project="masses")
    mass_weighted_Bx_map_ij = slice_gas(**common_arguments_xy, project="mass_weighted_Bx")
    mass_weighted_By_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_By")
    mass_weighted_Bz_map_ij  = slice_gas(**common_arguments_xy, project="mass_weighted_Bz")
    Bx_map_ij  = mass_weighted_Bx_map_ij /mass_map_ij 
    By_map_ij  = mass_weighted_By_map_ij /mass_map_ij 
    Bz_map_ij  = mass_weighted_Bz_map_ij /mass_map_ij 


    # Invert y axis to convert ij map to xy map
    Bx_map_xy = np.flipud(Bx_map_ij)
    By_map_xy = np.flipud(By_map_ij)
    Bz_map_xy = np.flipud(Bz_map_ij)

    B = np.stack((Bx_map_xy, By_map_xy, Bz_map_xy), axis=-1)
    return {'B':B, 'z':height}


def make_SPH_grid_B(data,
                    Nz = 2*res
                    ):
    """
    Rendering B field from SPH particle disrtibution
    Parameters: data: SPH snapshot data
                Nz: number of samples along z
    Returns: 
                Bx,By,Bz: 3 3d arrays of shape (N_z, N_x, N_y)
    """
    
    # Retrieve particle attributes of interest
    B = data.gas.magnetic_flux_densities
    #h = data.gas.smoothing_lengths
    #Np = len(h)

    # Generate mass weighted maps of quantities of interest
    data.gas.mass_weighted_Bx = data.gas.masses * B[:,0]
    data.gas.mass_weighted_By = data.gas.masses * B[:,1]
    data.gas.mass_weighted_Bz = data.gas.masses * B[:,2]

    slices = np.linspace(0,Lbox[2].value, Nz, endpoint=False)*Lbox[2].units
    #slices+=Lbox[2]/(2*Nz)

    Bx = []
    By = []
    Bz = []

    for slice in slices:
        Bx.append(make_2D_slice_z(data, slice)['B'][:,:,0])
        By.append(make_2D_slice_z(data, slice)['B'][:,:,1])
        Bz.append(make_2D_slice_z(data, slice)['B'][:,:,2])

    return np.array(Bx),np.array(By),np.array(Bz)


def compute_fourier_transform_1d(Q, dx):

    """
    Applying fast Fourier transformation to a quantity Q
    Parameters: Q: 1D array
                dx: spacial discretization
    Returns: 
                k: wavenumbers for each mode
                Q_k: Fourier image of Q
    """

    # Grid size
    N = len(Q)

    # Compute Fourier transforms
    #Q_k = np.fft.fftn(Q) * dx / Lbox[2] 
    Q_k = np.fft.fft(Q) * dx / Lbox[2] 

    # Compute the corresponding wavenumbers
    k = np.fft.fftfreq(N, d=dx) * 2 * np.pi
    
    return k, Q_k


def transform_B_along_z(Bx,
                        By,
                        Bz,
                        Nz=2*res
                        ):
    """
    Applying Fourier transformation along z to the B field array
    Parameters: Bx,By,Bz: 3D arrays with B field components
                Nz: number of samples along z
    Returns: 
                k: wavenumbers for each mode
                b_xy: array of shape (N_x,N_y, N_modes, N_dim=3) - Fourier image of B
    """

    dx = Lbox[2]/Nz

    b_xy = []
    for i in range(res):
        b_xy_row = []
        for j in range(res):
            k, bx_k = compute_fourier_transform_1d(Bx[:,i,j], dx)
            k, by_k = compute_fourier_transform_1d(By[:,i,j], dx)
            k, bz_k = compute_fourier_transform_1d(Bz[:,i,j], dx)
            b_k = np.stack((bx_k, by_k, bz_k), axis=-1)
            b_xy_row.append(b_k)
        b_xy.append(b_xy_row)

    return k, np.array(b_xy)


def save_to_sav(savedat:dict,
                outputfilename:str="b0perp.sav"
                ):
    """
    Save the data from dict to .sav file
    Parameters: savedat: data to save, dict of the form {'item1':array, ...}
                outputfilename: output filename
    Returns: 
    """

    # Save in IDL .sav format with the variable name 'b0perp'
    sp.io.savemat(outputfilename, savedat)


def load_cmap(img_addr:str,
              cmapname:str
              ):
    """
    Load custom colormaps (loads an image and creates colormap).
    Parameters: img_addr: adress of the colormap image
                cmapname: name of the colormap
    Returns: colormap object
    """

    img = imread(img_addr)

    colors_from_img = img[0, :, :]

    size = len(colors_from_img)

    res = LinearSegmentedColormap.from_list(cmapname, colors_from_img, N=size)
    return res


def generate_test_field(data
                        ):
    """
    Generate a test field for testing purposes.
    Parameters: data: snapshot data
    Returns: test_field: a 3D array
    """

    import numpy as np
    #res = 64 # 16 #64 #16
    Nz = 2*res

    # Define box size (to mimic Lbox from snapshot metadata)
    Lbox = data.metadata.boxsize
    
    #np.array([1.0, 1.0, 1.0])  # in arbitrary units

    # Create coordinates
    x = np.linspace(0, Lbox[0].value, res, endpoint=False) * Lbox.units
    y = np.linspace(0, Lbox[1].value, res, endpoint=False) * Lbox.units
    z = np.linspace(0, Lbox[2].value, Nz, endpoint=False) * Lbox.units

    Z, X, Y = np.meshgrid(z, x, y, indexing="ij")


    # Example magnetic fields
    Bx = np.zeros_like(Y)                          # uniform field along x
    By = np.zeros_like(Y)                               #np.zeros_like(Y)       #1.0*np.cos(2*np.pi*Z/Lbox[2])               # no y-component
    Bz = np.cos(2*np.pi*Z/Lbox[2]+X.value) #np.zeros_like(Z)                              #np.cos(2*np.pi*Z/Lbox[2])*np.cos(2*np.pi*X/Lbox[2]) #np.zeros_like(Y) #np.cos(2*np.pi*Z/Lbox[2]+0.5*np.pi*np.cos(2*np.pi*X/Lbox[0]))  #np.zeros_like(Z) #1.0*np.cos(2*np.pi*Z/Lbox[2])               # np.zeros_like(Z) 
    B_data = (Bx, By, Bz) * data.gas.magnetic_flux_densities.units
    print(Bx)
    return B_data


# Plotting

# Pick specific mode, nz
modez = 1 

# convert SPH data to grid
B_data = make_SPH_grid_B(data)

# Test plotter
#B_data = generate_test_field(data)

# Normalize B
B2 = B_data[0][:,:,:]**2+B_data[1][:,:,:]**2+B_data[2][:,:,:]**2
Brms = np.sqrt(np.mean(B2.flatten()))
print('B_rms = ',Brms)
B_data/=Brms

# Extract modes
k, b = transform_B_along_z(B_data[0],B_data[1],B_data[2], Nz=2*res)

# Print mode momentum
print('Mode: k=',k[modez])

# Calculate |b_i|, <b_x>,  <b_y>
b_abs = np.abs(b[:,:,modez,:])
bx_mean = np.mean(b[:,:,modez,0])
by_mean = np.mean(b[:,:,modez,1])
bx_phase = bx_mean/np.abs(bx_mean)
by_phase = by_mean/np.abs(by_mean)

# Re(<b_x>)>0
if np.real(bx_phase)<0:
    bx_phase*=-1
# Re(<b_y>)>0
if np.real(by_phase)<0:
    by_phase*=-1

# Display phase
print('Complex phase arg(b_x):',np.angle(bx_phase))
print('Complex phase arg(b_y):',np.angle(by_phase))

# Im(<b_x>)=0
theta = np.angle(b[:,:,modez,:]/bx_phase)
# Im(<b_y>)=0
theta_tilde = np.angle(b[:,:,modez,:]/by_phase)

# Fix phases to the interval \phi \in [-\pi,\pi]
theta[theta>np.pi] = np.pi - theta[theta>np.pi]
theta[theta<-np.pi] = theta[theta<-np.pi] - np.pi
theta_tilde[theta_tilde>np.pi] = np.pi - theta_tilde[theta_tilde>np.pi]
theta_tilde[theta_tilde<-np.pi] = theta_tilde[theta_tilde<-np.pi] - np.pi


# Generate figure
nx = 3
ny = 3
fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5 * nx, 5 * ny))

# Number of color levels
nlev=100 #17

# Colormaps
# Custom colormaps
#c1r = load_cmap('./cmap/c1r.png','c1r')
#c2 = load_cmap('./cmap/c2.png','c2')
# Select colormaps for r and theta
cmapr = "plasma" #c1r #"plasma"
cmaptheta = "twilight_shifted" #c1r #"RdBu_r" # "twilight_shifted"

# Add plots
a00 = ax[0, 0].contourf(
    periodic_continue(b_abs[:,:,0].T),
    cmap=cmapr,
    extend="both",
    levels=np.linspace(0, 1, nlev),
)
a01 = ax[0, 1].contourf(
    periodic_continue(b_abs[:,:,1].T),
    cmap=cmapr,
    extend="both",
    levels=np.linspace(0, 1, nlev),
)
a02 = ax[0, 2].contourf(
    periodic_continue(b_abs[:,:,2].T),
    cmap=cmapr,
    extend="both",
    levels=np.linspace(0, 1, nlev),
)
a10 = ax[1, 0].contourf(
    periodic_continue(theta[:,:,0].T),
    cmap=cmaptheta,
    extend="both",
    levels=np.linspace(-np.pi, np.pi, nlev),
)
a11 = ax[1, 1].contourf(
    periodic_continue(theta[:,:,1].T),
    cmap=cmaptheta,
    extend="both",
    levels=np.linspace(-np.pi, np.pi, nlev),
)
a12 = ax[1, 2].contourf(
    periodic_continue(theta[:,:,2].T),
    cmap=cmaptheta,
    extend="both",
    levels=np.linspace(-np.pi, np.pi, nlev),
)
a20 = ax[2, 0].contourf(
    periodic_continue(theta_tilde[:,:,0].T),
    cmap=cmaptheta,
    extend="both",
    levels=np.linspace(-np.pi, np.pi, nlev),
)
a21 = ax[2, 1].contourf(
    periodic_continue(theta_tilde[:,:,1].T),
    cmap=cmaptheta,
    extend="both",
    levels=np.linspace(-np.pi, np.pi, nlev),
)
a22 = ax[2, 2].contourf(
    periodic_continue(theta_tilde[:,:,2].T),
    cmap=cmaptheta,
    extend="both",
    levels=np.linspace(-np.pi, np.pi, nlev),
)

# Add ticks and labels
locs = [res/4 , res/2, 3 * res/4, res]
labels = [r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$',r'$2\pi$']
for ii in range(ny):
    ax[ii, 0].set_ylabel(r"$y$")
    ax[ii, 0].set_yticks(locs, labels)
for ii in range(ny):
    for jj in range(nx):
        ax[ii, jj].set_xlabel(r"$x$")
        ax[ii, jj].set_xticks(locs, labels)
        #ax[ii, jj].set_aspect("equal")
        ax[ii, jj].set_aspect("equal", adjustable="box")
        ax[ii, jj].set_xlim(0, res)
        ax[ii, jj].set_ylim(0, res)

# Add colorbar
ticks_r = [0,1]
ticks_phase = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
cbar1 = plt.colorbar(a00, ax=ax[0, 0], fraction=0.046, pad=0.04, ticks=ticks_r)
cbar1.set_label(r"$r_1(x,y) \;$")
cbar2 = plt.colorbar(a01, ax=ax[0, 1], fraction=0.046, pad=0.04, ticks=ticks_r)
cbar2.set_label(r"$r_2(x,y) \;$")
cbar3 = plt.colorbar(a02, ax=ax[0, 2], fraction=0.046, pad=0.04, ticks=ticks_r)
cbar3.set_label(r"$r_3(x,y) \;$")
cbar4 = plt.colorbar(a10, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=ticks_phase)
cbar4.set_label(r"$\theta_1(x,y) \;$")
cbar5 = plt.colorbar(a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=ticks_phase)
cbar5.set_label(r"$\theta_2(x,y) \;$")
cbar6 = plt.colorbar(a12, ax=ax[1, 2], fraction=0.046, pad=0.04, ticks=ticks_phase)
cbar6.set_label(r"$\theta_3(x,y) \;$")
cbar7 = plt.colorbar(a20, ax=ax[2, 0], fraction=0.046, pad=0.04, ticks=ticks_phase)
cbar7.set_label(r"$\tilde{\theta}_1(x,y) \;$")
cbar8 = plt.colorbar(a21, ax=ax[2, 1], fraction=0.046, pad=0.04, ticks=ticks_phase)
cbar8.set_label(r"$\tilde{\theta}_2(x,y) \;$")
cbar9 = plt.colorbar(a22, ax=ax[2, 2], fraction=0.046, pad=0.04, ticks=ticks_phase)
cbar9.set_label(r"$\tilde{\theta}_3(x,y) \;$")

fig.tight_layout()

plt.savefig(args.output, dpi=220)
