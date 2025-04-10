import numpy as np
import h5py
import unyt
import matplotlib.pyplot as plt
import argparse

from swiftsimio import load
from swiftsimio.visualisation.volume_render import render_gas


# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
argparser.add_argument("resolution")
args = argparser.parse_args()

# Load snapshot
filename = args.input
data = load(filename)

time = data.metadata.time
print("Plotting at %" % round(time,2))

Lbox = data.metadata.boxsize

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

# Retrieve particle attributes of interest
v = data.gas.velocities
B = data.gas.magnetic_flux_densities
h = data.gas.smoothing_lengths
minh = np.min(h.value)

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_vx = data.gas.masses * v[:,0]
data.gas.mass_weighted_vy = data.gas.masses * v[:,1]
data.gas.mass_weighted_vz = data.gas.masses * v[:,2]

data.gas.mass_weighted_Bx = data.gas.masses * B[:,0]
data.gas.mass_weighted_By = data.gas.masses * B[:,1]
data.gas.mass_weighted_Bz = data.gas.masses * B[:,2]

res = int(args.resolution)

common_arguments = dict(
    data=data, resolution=res, parallel=True,periodic=True,
)


mass_cube = render_gas(**common_arguments, project="masses")

mass_weighted_vx_cube = render_gas(**common_arguments, project="mass_weighted_vx")
mass_weighted_vy_cube = render_gas(**common_arguments, project="mass_weighted_vy")
mass_weighted_vz_cube = render_gas(**common_arguments, project="mass_weighted_vz")

mass_weighted_Bx_cube = render_gas(**common_arguments, project="mass_weighted_Bx")
mass_weighted_By_cube = render_gas(**common_arguments, project="mass_weighted_By")
mass_weighted_Bz_cube = render_gas(**common_arguments, project="mass_weighted_Bz")

vx_cube = mass_weighted_vx_cube/mass_cube
vy_cube = mass_weighted_vy_cube/mass_cube
vz_cube = mass_weighted_vz_cube/mass_cube

Bx_cube = mass_weighted_Bx_cube/mass_cube
By_cube = mass_weighted_By_cube/mass_cube
Bz_cube = mass_weighted_Bz_cube/mass_cube

def compute_magnetic_power_spectrum(Qx, Qy, Qz, dx,nbins=res):
    # Grid size
    Nx, Ny, Nz = Qx.shape

    # Compute Fourier transforms
    Qx_k = np.fft.fftn(Qx)
    Qy_k = np.fft.fftn(Qy)
    Qz_k = np.fft.fftn(Qz)

    # Compute power spectrum (squared magnitude of Fourier components)
    Q_power_k = np.abs(Qx_k)**2 + np.abs(Qy_k)**2 + np.abs(Qz_k)**2

    # Compute the corresponding wavenumbers
    kx = np.fft.fftfreq(Nx, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(Ny, d=dx) * 2 * np.pi
    kz = np.fft.fftfreq(Nz, d=dx) * 2 * np.pi

    # Create 3D arrays of wavevectors
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Flatten arrays for binning
    k_mag_flat = k_mag.flatten()
    Q_power_flat = Q_power_k.flatten()

    # Define k bins (you can tweak bin size)
    k_bins = np.linspace(0, np.max(k_mag), num=nbins)
    k_bin_centers = 0.5 * (k_bins[1:] + k_bins[:-1])

    # Bin the power spectrum
    power_spectrum, _ = np.histogram(k_mag_flat, bins=k_bins, weights=Q_power_flat)
    counts, _ = np.histogram(k_mag_flat, bins=k_bins)
    
    # Avoid division by zero
    power_spectrum = np.where(counts > 0, power_spectrum / counts, 0)

    return k_bin_centers, power_spectrum

# plot magnetic field spectrum
ks, Pb = compute_magnetic_power_spectrum(Bx_cube.value,By_cube.value,Bz_cube.value, dx = Lbox[0].value/(res) )
# plot velocity spectrum
ks, Pv = compute_magnetic_power_spectrum(vx_cube.value,vy_cube.value,vz_cube.value, dx = Lbox[0].value/(res) )


plt.plot(ks,Pb,color='red',linestyle='dashed',label='$P_B(k)$')
plt.plot(ks,Pv,color='blue',label='$P_v(k)$')

# plot spectral lines
ksmock = np.logspace(1,2,10)
p = -5/3
plt.plot(ksmock,1e10*ksmock**(p),color='black',linestyle='dashed',label = '$k^{-5/3}$')

# plot spectral lines
ksmock = np.logspace(-0.5,0.5,10)
p = 3/2
plt.plot(ksmock,1e6*ksmock**(p),color='black',linestyle='dashed',label = '$k^{3/2}$')

plt.axvline(x=2*np.pi/minh, color='gray',linestyle='dashed',label = '$k_{res}$')

plt.yscale('log')
plt.xscale('log')
plt.xlabel('k')
plt.ylabel('P(k)')
plt.grid()
plt.legend()
#plt.xlim([0,15])
plt.savefig(args.output)




