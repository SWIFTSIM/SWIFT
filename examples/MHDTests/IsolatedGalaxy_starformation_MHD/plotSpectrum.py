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
#argparser.add_argument("resolution")
args = argparser.parse_args()

# Load snapshot
filename = args.input
data = load(filename)

time = np.round(data.metadata.time,2)
print("Plotting at ", time)

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
rho = data.gas.densities
B = data.gas.magnetic_flux_densities
h = data.gas.smoothing_lengths
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)

divB = data.gas.magnetic_divergences
minh = np.min(h.value)

# Renormalize quantities
mu0 = 1.25663706127e-1 * unyt.g * unyt.cm / (unyt.s ** 2 * unyt.statA ** 2)
v = v * (rho[:, None]/2)**0.5
B = B / np.sqrt(2*mu0)

#Npside = int(len(h)**(1/3))+1

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_vx = data.gas.masses * v[:,0]
data.gas.mass_weighted_vy = data.gas.masses * v[:,1]
data.gas.mass_weighted_vz = data.gas.masses * v[:,2]

data.gas.mass_weighted_Bx = data.gas.masses * B[:,0]
data.gas.mass_weighted_By = data.gas.masses * B[:,1]
data.gas.mass_weighted_Bz = data.gas.masses * B[:,2]

data.gas.mass_weighted_error = data.gas.masses * np.maximum(h * abs(divB) / (normB + 0.01 * np.max(normB)), 1e-6)

#res = 128 #Npside #int(args.resolution)

#Lslice = 0.5*Lbox
Lslice_kPc = 20  # 2*21.5
Lslice = Lslice_kPc * 3.08e18 * 1e3 * unyt.cm

k_slice = 2*np.pi/Lslice
k_res = np.max(2*np.pi/h)

res = int((2*Lslice/np.min(h)).value) #Npside #int(args.resolution)

# resolution limiter
resmax=512
print('Suggested maximal resoluiton: ',res)
res = min([res,resmax])

print('Spectral resolution: ',res)

center = Lbox/2
visualise_region = [
    center[0] - Lslice,
    center[0] + Lslice,
    center[1] - Lslice,
    center[1] + Lslice,
    center[2] - Lslice,
    center[2] + Lslice,
]

common_arguments = dict(
    data=data, resolution=res, parallel=True,periodic=True,region=visualise_region,
)


mass_cube = render_gas(**common_arguments, project="masses")
min_nonzero_m =np.min( mass_cube.flatten()[mass_cube.flatten().value!=0.0])
mass_cube[mass_cube.value==0]=1e-2*min_nonzero_m

mass_weighted_vx_cube = render_gas(**common_arguments, project="mass_weighted_vx")
mass_weighted_vy_cube = render_gas(**common_arguments, project="mass_weighted_vy")
mass_weighted_vz_cube = render_gas(**common_arguments, project="mass_weighted_vz")

mass_weighted_Bx_cube = render_gas(**common_arguments, project="mass_weighted_Bx")
mass_weighted_By_cube = render_gas(**common_arguments, project="mass_weighted_By")
mass_weighted_Bz_cube = render_gas(**common_arguments, project="mass_weighted_Bz")

mass_weighted_error_cube = render_gas(**common_arguments, project="mass_weighted_error")


vx_cube = mass_weighted_vx_cube/mass_cube
vy_cube = mass_weighted_vy_cube/mass_cube
vz_cube = mass_weighted_vz_cube/mass_cube

Bx_cube = mass_weighted_Bx_cube/mass_cube
By_cube = mass_weighted_By_cube/mass_cube
Bz_cube = mass_weighted_Bz_cube/mass_cube

error_cube = mass_weighted_error_cube/mass_cube 

unit_energy = unyt.g * unyt.cm**2 / unyt.s**2
unit_energy_density = unit_energy / unyt.cm**3
unit_sqrt_energy_density = (unit_energy_density)**0.5
unit_length = 1e3*unyt.pc

vx_cube.convert_to_units(unit_sqrt_energy_density)
vy_cube.convert_to_units(unit_sqrt_energy_density)
vz_cube.convert_to_units(unit_sqrt_energy_density)

Bx_cube.convert_to_units(unit_sqrt_energy_density)
By_cube.convert_to_units(unit_sqrt_energy_density)
Bz_cube.convert_to_units(unit_sqrt_energy_density)

k_slice = 2*np.pi/Lslice
k_res = np.max(2*np.pi/h)

k_slice.convert_to_units(1/unit_length)
k_res.convert_to_units(1/unit_length)

def compute_power_spectrum_scal(Q, dx, nbins):
    # Grid size
    Nx, Ny, Nz = Q.shape

    # Compute Fourier transforms
    Q_k = np.fft.fftn(Q)

    # Compute power spectrum (squared magnitude of Fourier components)
    Q_power_k = np.abs(Q_k)**2

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
    k_bins = np.linspace(0, np.max(2*np.pi/minh), num=nbins)

    # exclude 0th mode:
    k_bins = k_bins[1:]
    k_mag_flat = k_mag_flat[1:]
    Q_power_flat = Q_power_flat[1:]

    k_bin_centers = 0.5 * (k_bins[1:] + k_bins[:-1])

    # Bin the power spectrum
    power_spectrum, _ = np.histogram(k_mag_flat, bins=k_bins, weights=Q_power_flat)
    counts, _ = np.histogram(k_mag_flat, bins=k_bins)
    
    # Avoid division by zero
    power_spectrum = np.where(counts > 0, power_spectrum / counts, 0)

    return k_bin_centers, power_spectrum

def compute_power_spectrum_vec(Qx, Qy, Qz, dx,nbins):
    # Ensure all arrays are unyt arrays with same units
    dx = dx #.to(unyt.pc)
    volume_element = dx**3  # single cell volume in cm^3

    # Get shape
    Nx, Ny, Nz = Qx.shape

    # Fourier transform (keep units)
    Qx_k = np.fft.fftn(Qx.value) * volume_element.value
    Qy_k = np.fft.fftn(Qy.value) * volume_element.value
    Qz_k = np.fft.fftn(Qz.value) * volume_element.value

    # Reapply units
    Qx_k = Qx_k * Qx.units * volume_element.units
    Qy_k = Qy_k * Qy.units * volume_element.units
    Qz_k = Qz_k * Qz.units * volume_element.units

    # Power in k-space
    Pk = (np.abs(Qx_k)**2 + np.abs(Qy_k)**2 + np.abs(Qz_k)**2)#.to(unyt.g/(unyt.s**2)*unyt.cm**5)  # or appropriate units

    # Compute wavenumbers
    kx = np.fft.fftfreq(Nx, d=dx.value) * 2 * np.pi / dx.units #unyt.cm**-1
    ky = np.fft.fftfreq(Ny, d=dx.value) * 2 * np.pi / dx.units #unyt.cm**-1
    kz = np.fft.fftfreq(Nz, d=dx.value) * 2 * np.pi / dx.units #unyt.cm**-1
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Flatten and bin
    k_flat = k_mag.flatten() #.to('1/cm')
    Pk_flat = Pk.flatten() #.to(unyt.g/(unyt.s**2)*unyt.cm**5)  # this should be correct after unit reduction

    # Bin in k-space
    k_min_log = np.log10(np.min(k_flat.value[k_flat.value!=0]))
    k_max_log = np.log10(np.max(k_flat.value))
    k_bins = np.logspace(k_min_log,k_max_log, num=nbins) * k_flat.units  #*unyt.cm**-1

    # exclude 0th mode:
    k_bins = k_bins[1:]
    k_flat = k_flat[1:]
    Pk_flat = Pk_flat[1:]

    # converting to Energy per unit wavevector
    Ek_flat = 4*np.pi * k_flat**2 * Pk_flat

    k_bin_centers = 0.5 * (k_bins[1:] + k_bins[:-1])

    power_spectrum = np.zeros(len(k_bin_centers)) * Ek_flat.units
    counts = np.zeros(len(k_bin_centers))
   
    for i in range(len(k_bin_centers)):
        in_bin = (k_flat >= k_bins[i]) & (k_flat < k_bins[i+1])
        power_spectrum[i] = Ek_flat[in_bin].sum()
        counts[i] = in_bin.sum()

    # Normalize per mode (optional)
    power_spectrum = np.where(counts > 0, power_spectrum / counts, 0 * power_spectrum.units)

    # mask non-zero
    mask_zeros = power_spectrum==0
    k_bin_centers = k_bin_centers[~mask_zeros]
    power_spectrum = power_spectrum[~mask_zeros]

    return k_bin_centers, power_spectrum

dx = 2*Lslice.to(unit_length)/(res)

# plot magnetic field spectrum
ks, Eb = compute_power_spectrum_vec(Bx_cube,By_cube,Bz_cube, dx = dx, nbins=res-1 )
# plot velocity spectrum
ks, Ev = compute_power_spectrum_vec(vx_cube,vy_cube,vz_cube, dx = dx, nbins=res-1 )
# plot divergence error spectrum
#ks, Perr = compute_power_spectrum_scal(error_cube, dx = dx, nbins=res-1 )

Eb.convert_to_units(unyt.erg*(1e3*unyt.pc))
Ev.convert_to_units(unyt.erg*(1e3*unyt.pc))
ks.convert_to_units(1e-3/unyt.pc)


fig, ax = plt.subplots(figsize=(10, 6.2))

ax.plot(ks,Ev.value,color='blue',label='$E_v(k)$')
ax.plot(ks,Eb.value,color='red',linestyle='solid',label='$E_B(k)$')
#ax.plot(ks,Perr,color='purple',label='$P_{R_{0}}(k)$')

# resolution line
ax.axvline(x=k_res, color='brown',linestyle='solid',label = r'$k_{\mathrm{res}}$')

# self gravity softening line
ksoft = 2*np.pi / (0.2 * 1e3 * unyt.pc)
ksoft.convert_to_units(1e-3/unyt.pc)
ax.axvline(x=ksoft, color='brown',linestyle='dashdot',label = r'$k_{\mathrm{soft}}$')

# plot spectral lines
ksmock = np.logspace(np.log10(ksoft.value)-1,np.log10(ksoft.value)-0.5,10)*ksoft.units
p = -5/3
Eright = 1e57 * unyt.erg*(1e3*unyt.pc)
kright = ksoft
kright.convert_to_units(1e-3/unyt.pc)
Emock = Eright*(ksmock/kright)**(p)
Emock.convert_to_units(unyt.erg*(1e3*unyt.pc))
ax.plot(ksmock,Emock,color='black',linestyle='dashed')#,label = '$k^{-5/3}$')

klabel = ksoft/10**0.45
klabel.convert_to_units(1e-3/unyt.pc)
Elabel = Eright*(klabel/kright)**(p)
Elabel.convert_to_units(unyt.erg*(1e3*unyt.pc))
ax.text(klabel,Elabel,r'$k^{-5/3}$')

# plot spectral lines
ksmock = np.logspace(np.log10(k_slice.value),np.log10(k_slice.value)+0.5,10)*k_slice.units
p = 3/2
Eleft = 1e57 * unyt.erg*(1e3*unyt.pc)
kleft = k_slice
kleft.convert_to_units(1e-3/unyt.pc)
Emock = Eleft*(ksmock/kleft)**(p)
Emock.convert_to_units(unyt.erg*(1e3*unyt.pc))
ax.plot(ksmock,Emock,color='black',linestyle='dashed')#,label = '$k^{-5/3}$')

klabel = k_slice/10**0.15
klabel.convert_to_units(1e-3/unyt.pc)
Elabel = Eleft*(klabel/kleft)**(p)
Elabel.convert_to_units(unyt.erg*(1e3*unyt.pc))
ax.text(klabel,Elabel,r'$k^{3/2}$')

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$\mathrm{k}$ $[\mathrm{kpc}^{-1}]$', fontsize=30)
ax.set_ylabel(r'$\mathrm{P}(\mathrm{k})$ $[\mathrm{erg}\cdot\mathrm{kpc}]$', fontsize=30)
ax.tick_params(axis='both', labelsize=20)
#ax.grid()
ax.legend()
fig.tight_layout()
plt.savefig(args.output)

