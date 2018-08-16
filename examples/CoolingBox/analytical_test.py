import sys
import matplotlib
matplotlib.use("Agg")
from pylab import *
import h5py

# Plot parameters
params = {'axes.labelsize': 10,
'axes.titlesize': 10,
'font.size': 12,
'legend.fontsize': 12,
'xtick.labelsize': 10,
'ytick.labelsize': 10,
'text.usetex': True,
 'figure.figsize' : (3.15,3.15),
'figure.subplot.left'    : 0.2,
'figure.subplot.right'   : 0.99,
'figure.subplot.bottom'  : 0.2,
'figure.subplot.top'     : 0.9,
'figure.subplot.wspace'  : 0.15,
'figure.subplot.hspace'  : 0.12,
'lines.markersize' : 6,
'lines.linewidth' : 3.,
'text.latex.unicode': True
}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Times']})


import numpy as np
import h5py as h5
import sys

# function for interpolating 1d table of cooling rates
def interpol_lambda(u_list,cooling_rate,u):
	if u < u_list[0]: 
		return[cooling_rate[0]]
	if u > u_list[len(u_list)-1]: 
		return[cooling_rate[len(cooling_rate)-1]]
	j = 0
	while u_list[j+1] < u:
		j += 1
	cooling = cooling_rate[j]*(u_list[j+1]-u)/(u_list[j+1]-u_list[j]) + cooling_rate[j+1]*(u-u_list[j])/(u_list[j+1]-u_list[j])
	return cooling

# File containing the total energy
stats_filename = "./energy.txt"

# First snapshot
snap_filename = "coolingBox_0000.hdf5"

# Some constants in cgs units
k_b = 1.38E-16 #boltzmann
m_p = 1.67e-24 #proton mass

# Initial conditions set in makeIC.py
T_init = 1.0e7

# Read the initial state of the gas
f = h5.File(snap_filename,'r')
rho = np.mean(f["/PartType0/Density"])
pressure = np.mean(f["/PartType0/Pressure"])

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]
unit_vel = unit_length/unit_time

rho = rho*unit_mass/(unit_length**3)

# Read snapshots
if len(sys.argv) >= 4:
	nsnap = int(sys.argv[5])
else:
	print("Need to specify number of snapshots, defaulting to 100.")
	nsnap = 100
npart = 4096
u_snapshots_cgs = zeros(nsnap)
u_part_snapshots_cgs = zeros((nsnap,npart))
t_snapshots_cgs = zeros(nsnap)
scale_factor = zeros(nsnap)
for i in range(nsnap):
    snap = h5.File("coolingBox_%0.4d.hdf5"%i,'r')
    u_part_snapshots_cgs[i][:] = snap["/PartType0/InternalEnergy"][:]*(unit_length**2)/(unit_time**2)
    t_snapshots_cgs[i] = snap["/Header"].attrs["Time"] * unit_time
    scale_factor[i] = snap["/Header"].attrs["Scale-factor"]

# calculate analytic solution. Since cooling rate is constant,
# temperature decrease in linear in time.
# read Lambda and temperatures from table
internal_energy = []
cooling_rate = []
file_in = open('cooling_output.dat', 'r')
for line in file_in:
        data = line.split()
        internal_energy.append(float(data[0]))
        cooling_rate.append(-float(data[1]))

tfinal = t_snapshots_cgs[nsnap-1]
tfirst = t_snapshots_cgs[0]
nt = nsnap*10
dt = (tfinal-tfirst)/nt

t_sol = np.zeros(int(nt))
u_sol = np.zeros(int(nt))
lambda_sol = np.zeros(int(nt))
u = np.mean(u_part_snapshots_cgs[0,:])/scale_factor[0]**2
u_sol[0] = u
t_sol[0] = tfirst
# integrate explicit ODE
for j in range(int(nt-1)):
	t_sol[j+1] = t_sol[j] + dt
	Lambda_net = interpol_lambda(internal_energy,cooling_rate,u_sol[j])
	if j < 10:
		print(u,Lambda_net)
	if int(sys.argv[4]) == 1:
		nH = 0.702*rho/(m_p)/scale_factor[0]**3
		ratefact = nH*0.702/m_p
	else:
		nH = 0.752*rho/(m_p)/scale_factor[0]**3
		ratefact = nH*0.752/m_p
	u_next = u - Lambda_net*ratefact*dt
	u_sol[j+1] = u_next
	u = u_next
u_sol = u_sol*scale_factor[0]**2

# swift solution
mean_u = np.zeros(nsnap)
for j in range(nsnap):
	mean_u[j] = np.mean(u_part_snapshots_cgs[j,:])

# plot and save results
log_nh = float(sys.argv[2])
redshift = float(sys.argv[1])
p = plt.figure(figsize = (6,6))
p1, = plt.loglog(t_sol,u_sol,linewidth = 0.5,color = 'k', marker = '*', markersize = 1.5,label = 'explicit ODE')
p2, = plt.loglog(t_snapshots_cgs,mean_u,linewidth = 0.5,color = 'r', marker = '*', markersize = 1.5,label = 'swift mean alexei')
l = legend(handles = [p1,p2])
xlabel("${\\rm{Time~[s]}}$", labelpad=0, fontsize = 14)
ylabel("Internal energy ${\\rm{[erg \cdot g^{-1}]}}$", fontsize = 14)
if int(sys.argv[4]) == 1:
	title('$n_h = 10^{' + "{}".format(log_nh) + '} \\rm{cm}^{-3}, z = ' + "{}".format(redshift) + '$, solar metallicity,\n relative error alexei: ' + "{0:.4f}".format( (u_sol[int(nt)-1] - mean_u[nsnap-1])/(u_sol[int(nt)-1])), fontsize = 16)
	name = "z_"+str(sys.argv[1])+"_nh_"+str(sys.argv[2])+"_pressure_"+str(sys.argv[3])+"_solar.png"
elif int(sys.argv[4]) == 0:
	title('$n_h = 10^{' + "{}".format(log_nh) + '} \\rm{cm}^{-3}, z = ' + "{}".format(redshift) + '$, zero metallicity,\n relative error alexei: ' + "{0:.4f}".format( (u_sol[int(nt)-1] - mean_u[nsnap-1])/(u_sol[int(nt)-1])), fontsize = 16)
	name = "z_"+str(sys.argv[1])+"_nh_"+str(sys.argv[2])+"_pressure_"+str(sys.argv[3])+"_zero_metal.png"

savefig(name, dpi=200)

print('Final internal energy ode, swift, error ' + str(u_sol[int(nt)-1]) + ' ' + str(mean_u[nsnap-1])  + ' ' + str( (u_sol[int(nt)-1] - mean_u[nsnap-1])/(u_sol[int(nt)-1])))
