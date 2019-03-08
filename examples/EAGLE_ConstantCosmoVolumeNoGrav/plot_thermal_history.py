import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py as h5
import swiftsimio
import sys
import glob
import numpy as np

## read command line arguments
snapshot_name = sys.argv[1]

params = {'axes.labelsize': 16,
'axes.titlesize': 16,
'font.size': 16,
'legend.fontsize': 12,
'xtick.labelsize': 16,
'ytick.labelsize': 16,
'text.usetex': False,
'figure.figsize' : (9,6),
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
schaye_T_lower_error = data_schaye[3] - data_schaye[4]
schaye_T_upper_error = data_schaye[5] - data_schaye[3]
walther_T_lower_error = data_walther[1] - data_walther[2]
walther_T_upper_error = data_walther[3] - data_walther[1]

############### Read in simulation data

## First, get list of all snapshots
reg_exp = "%s*.hdf5" %snapshot_name
snap_list = glob.glob(reg_exp)

z = np.array([])
T_mean = np.array([])
T_std = np.array([])

### loop through list
for snap in snap_list:
    
    # This loads all metadata but explicitly does _not_ read any particle data yet
    data = swiftsimio.load(snap)

    # Get the redshift
    z = np.append(z,data.metadata.z)
    
    # Convert gas temperatures to right units
    data.gas.temperature.convert_to_cgs()

    # Get mean and standard deviation of temperature
    T_mean = np.append(T_mean,np.mean(data.gas.temperature))
    T_std = np.append(T_std,np.std(data.gas.temperature))

T_mean /= 1.0e4
T_std /= 1.0e4

### sort arrays into redshift order
ind_sorted = np.argsort(z)
z = z[ind_sorted]
T_mean = T_mean[ind_sorted]
T_std = T_std[ind_sorted]

### from the first snapshot, get code information
d = swiftsimio.load(snap_list[0])
code_info = d.metadata.code
git_branch = code_info["Git Branch"].decode('UTF-8')
git_revision = code_info["Git Revision"].decode('UTF-8')
plot_title = git_branch + "    " + git_revision


# Make plot of temperature evolution  --------------------------------
fig, ax = plt.subplots()
ax.fill_between(z,T_mean - T_std,T_mean + T_std,alpha = 0.5)
ax.plot(z,T_mean,label = "Simulation")
ax.errorbar(data_schaye[0],data_schaye[3], xerr = [schaye_z_lower_error,schaye_z_upper_error],yerr = [schaye_T_lower_error,schaye_T_upper_error], fmt = 'ko', label = "Schaye+ 2000",zorder = 20,capsize = 4.0,capthick = 1.0,alpha = 0.9)
ax.errorbar(data_walther[0],data_walther[1],yerr = [walther_T_lower_error,walther_T_upper_error], fmt = 'rd', label = "Walther+ 2019",zorder = 30,capsize = 4.0,capthick = 1.0,alpha = 0.7)
ax.set_xlim(0.0,15.0)
ax.set_ylim(0.0,3.0)
ax.set_xlabel("Redshift")
ax.set_ylabel(r"$T\,/\,10^4K$")
ax.set_title(plot_title)
ax.legend(loc = 0)
fig.tight_layout()
fig.savefig("thermal_history.pdf",format = "pdf")
