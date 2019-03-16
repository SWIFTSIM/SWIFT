import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py as h5
import swiftsimio
import sys
import glob
import unyt
import numpy as np

## read command line arguments
snapshot_name = sys.argv[1]

params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 9,
'legend.fontsize': 9,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': False,
'figure.figsize' : (4.15,3.15),
'figure.subplot.left'    : 0.12,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.12,
'figure.subplot.top'     : 0.99,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 2.,
'text.latex.unicode': True
}
plt.rcParams.update(params)

################# Read in observational data

data_schaye = np.genfromtxt("./datasets/schaye_et_al_2000_thermal_history.dat",skip_header = 2)
data_walther = np.genfromtxt("./datasets/walther_et_al_2019_thermal_history.dat",skip_header = 2)

data_schaye = data_schaye.T
data_walther = data_walther.T

schaye_z_lower_error = data_schaye[0] - data_schaye[1]
schaye_z_upper_error = data_schaye[2] - data_schaye[0]
schaye_T_lower_error = np.log10(data_schaye[3]*1.0e4) - np.log10(data_schaye[4]*1.0e4)
schaye_T_upper_error = np.log10(data_schaye[5]*1.0e4) - np.log10(data_schaye[3]*1.0e4)
walther_T_lower_error = np.log10(data_walther[1]*1.0e4) - np.log10(data_walther[2]*1.0e4)
walther_T_upper_error = np.log10(data_walther[3]*1.0e4) - np.log10(data_walther[1]*1.0e4)

############### Read in simulation data

## First, get list of all snapshots
reg_exp = "%s*.hdf5" %snapshot_name
snap_list = glob.glob(reg_exp)

z = []
T_mean = []
T_std = []
rho_mean = []
rho_std = []

## loop through list
for snap in snap_list:
    
    # This loads all metadata but explicitly does _not_ read any particle data yet
    data = swiftsimio.load(snap)

    # Get the redshift
    z = np.append(z, data.metadata.z)
        
    # Convert gas temperatures and densities to right units
    data.gas.temperature.convert_to_cgs()

    # Get mean and standard deviation of temperature
    T_mean.append(np.mean(data.gas.temperature) * data.gas.temperature.units)
    T_std.append(np.std(data.gas.temperature) * data.gas.temperature.units)

    # Get mean and standard deviation of density
    rho_mean.append(np.mean(data.gas.density) * data.gas.density.units)
    rho_std.append(np.std(data.gas.density) * data.gas.density.units)
    
## Turn into numpy arrays
T_mean = np.array(T_mean) * data.gas.temperature.units
T_std = np.array(T_std) * data.gas.temperature.units
rho_mean = np.array(rho_mean) * data.gas.density.units
rho_std = np.array(rho_std) * data.gas.density.units

## Put Density into units of mean baryon density today

# first calculate rho_bar_0 from snapshot metadata
### from the first snapshot, get cosmology information
d = swiftsimio.load(snap_list[0])
cosmology = d.metadata.cosmology
H0 = cosmology["H0 [internal units]"] / (d.units.time)
Omega_bar = cosmology["Omega_b"]

### now calculate rho_bar_0 and divide through
rho_bar_0 = 3.0 * H0**2 / (8. * np.pi * unyt.G) * Omega_bar 
rho_mean /= rho_bar_0
rho_std /= rho_bar_0

### sort arrays into redshift order
ind_sorted = np.argsort(z)
z = z[ind_sorted]
T_mean = T_mean[ind_sorted]
T_std = T_std[ind_sorted]
rho_mean = rho_mean[ind_sorted]
rho_std = rho_std[ind_sorted]

### from the first snapshot, get code information
d = swiftsimio.load(snap_list[0])
code_info = d.metadata.code
git_branch = code_info["Git Branch"].decode('UTF-8')
git_revision = code_info["Git Revision"].decode('UTF-8')
hydro_metadata = d.metadata.hydro_scheme
scheme = hydro_metadata["Scheme"].decode('UTF-8')
params = d.metadata.parameters

subgrid_metadata = d.metadata.subgrid_scheme
cooling_model = subgrid_metadata["Cooling Model"].decode('UTF-8')
chemistry_model = subgrid_metadata["Chemistry Model"].decode('UTF-8')

if cooling_model == 'EAGLE':
    z_r_H = float(params['EAGLECooling:H_reion_z'])
    H_heat_input = float(params['EAGLECooling:H_reion_eV_p_H'])
    z_r_He_centre = float(params['EAGLECooling:He_reion_z_centre'])
    z_r_He_sigma = float(params['EAGLECooling:He_reion_z_sigma'])
    He_heat_input = float(params['EAGLECooling:He_reion_eV_p_H'])

metallicity = "Unknown"
if chemistry_model == 'EAGLE':
    metallicity = float(params['EAGLEChemistry:init_abundance_metal'])
    
# Make plot of temperature evolution  --------------------------------
fig = plt.figure()

# Plot sim properties
if cooling_model == 'EAGLE':
    plt.plot([z_r_H, z_r_H], [3.4, 4.4], 'k--', alpha=0.5, lw=0.7)
    plt.text(z_r_H + 0.1, 3.55, "H reion.", rotation=90, alpha=0.5, fontsize=7, va="bottom")
    plt.plot([z_r_He_centre, z_r_He_centre], [3.4, 4.4], 'k--', alpha=0.5, lw=0.7)
    plt.text(z_r_He_centre + 0.1, 3.55, "HeII reion.", rotation=90, alpha=0.5, fontsize=7, va="bottom")
    
# Plot observational data
plt.errorbar(data_schaye[0],
             np.log10(data_schaye[3]*1.0e4),
             xerr = [schaye_z_lower_error,schaye_z_upper_error],
             yerr = [schaye_T_lower_error,schaye_T_upper_error],
             fmt='s', mec='0.3', color='0.3', markersize=4, markeredgewidth=0.5, linewidth=0.5, mfc='w', label="Schaye et al. (2000)")
plt.errorbar(data_walther[0],
             np.log10(data_walther[1]*1.0e4),
             yerr = [walther_T_lower_error,walther_T_upper_error],
              fmt='.', mec='0.3', color='0.3', markersize=7, markeredgewidth=0.5, linewidth=0.5, mfc='w', label = "Walther et al. (2019)")

# Plot simulation
plt.plot(z, np.log10(T_mean))

# Legend
plt.legend(loc="upper right", frameon=True, fontsize=8, handletextpad=0.1, facecolor="w", edgecolor="w", framealpha=1.)
plt.text(0.2, 4.8, "SWIFT %s \nCooling model: %s \n$\\rho=%.3f$ $\\Omega_{b}\\rho_{crit,0}$\n$Z=%s$"%(git_revision, cooling_model, rho_mean[-1], metallicity), va="top", ha="left", fontsize=8)

plt.xlim(0, 12.2)
plt.ylim(3.5,4.85)
plt.xlabel("Redshift", labelpad=-0.5)
plt.ylabel(r"$\log_{10}(T/K)$", labelpad=0)
plt.savefig("Temperature_evolution.png", dpi=200)


# Make plot of denisty evolution  --------------------------------
plt.rcParams.update({'figure.subplot.left'    : 0.14})
fig = plt.figure()

plt.text(0.2, 1.011, "SWIFT %s"%git_revision, va="top", ha="left", fontsize=8)

plt.fill_between(z,rho_mean - rho_std,rho_mean + rho_std,alpha = 0.5)
plt.plot(z,rho_mean)

plt.axhline(y = 1.0, linestyle = '--', color='k', alpha=0.5, lw=1.)

plt.xlim(0.0,12.2)
plt.ylabel(r"$\delta_b = \rho / \Omega_b\rho_{crit,0}$", labelpad=0.)
plt.ylim(0.988,1.012)
plt.yticks([0.99, 1., 1.01])
plt.xlabel("Redshift", labelpad=-0.5)
plt.savefig("Density_evolution.png", dpi=200)
