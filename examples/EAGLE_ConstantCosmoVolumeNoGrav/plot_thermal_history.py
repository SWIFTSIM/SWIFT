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
plot_filename = sys.argv[2]

params = {'axes.labelsize': 16,
'axes.titlesize': 12,
'font.size': 16,
'legend.fontsize': 12,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': False,
          'figure.figsize' : (9,9),
'figure.subplot.left'    : 0.1,
'figure.subplot.right'   : 0.95  ,
'figure.subplot.bottom'  : 0.15  ,
'figure.subplot.top'     : 0.9  ,
'figure.subplot.wspace'  : 0.0  ,
'figure.subplot.hspace'  : 0.15  ,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
plt.rcParams.update(params)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Times']})

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
rho_std = []### loop through list
for snap in snap_list:
    
    # This loads all metadata but explicitly does _not_ read any particle data yet
    data = swiftsimio.load(snap)

    # Get the redshift
    z = np.append(z,data.metadata.z)
    
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
params = d.metadata.parameters
z_r_H = float(params['EAGLECooling:H_reion_z'])
H_heat_input = float(params['EAGLECooling:H_reion_eV_p_H'])
z_r_He_centre = float(params['EAGLECooling:He_reion_z_centre'])
z_r_He_sigma = float(params['EAGLECooling:He_reion_z_sigma'])
He_heat_input = float(params['EAGLECooling:He_reion_eV_p_H'])

version_info= git_branch +'/'+ git_revision + '\n'
reion_info = "$z_{r,H} = %1.1f \; \Delta u_H = %1.1f \; z_{r,He,mid} = %1.1f \; z_{r,He,\sigma} = %1.1f \; \Delta u_{He} = %1.1f \; eV/m_H$" %(z_r_H,H_heat_input,z_r_He_centre,z_r_He_sigma,He_heat_input)

plot_title = version_info + reion_info

# Make plot of temperature evolution  --------------------------------
fig, axes = plt.subplots(2,1,sharex = True)
axes[0].fill_between(z,np.log10(T_mean - T_std),np.log10(T_mean + T_std),alpha = 0.5)
axes[0].plot(z,np.log10(T_mean),label = "Simulation")
axes[0].errorbar(data_schaye[0],np.log10(data_schaye[3]*1.0e4), xerr = [schaye_z_lower_error,schaye_z_upper_error],yerr = [schaye_T_lower_error,schaye_T_upper_error], fmt = 'ko', label = "Schaye+ 2000",zorder = 20,capsize = 4.0,capthick = 1.0,alpha = 0.9)
axes[0].errorbar(data_walther[0],np.log10(data_walther[1]*1.0e4),yerr = [walther_T_lower_error,walther_T_upper_error], fmt = 'rd', label = "Walther+ 2019",zorder = 30,capsize = 4.0,capthick = 1.0,alpha = 0.7)
axes[1].fill_between(z,rho_mean - rho_std,rho_mean + rho_std,alpha = 0.5)
axes[1].plot(z,rho_mean,label = "Simulation")
axes[1].axhline(y = 1.0,linestyle = '--',color = 'r')
axes[1].set_xlim(0.0,15.0)
axes[0].set_ylim(2.0,5.0)
axes[1].set_ylim(0.99,1.01)
axes[1].set_xlabel("Redshift")
axes[0].set_ylabel(r"$\log_{10}(T/K)$")
axes[1].set_ylabel(r"$\delta_b$")
axes[0].set_title(plot_title)
axes[0].legend(loc = 0)
axes[1].legend(loc = 0)
fig.tight_layout()
fig.savefig(plot_filename,format = "pdf")
