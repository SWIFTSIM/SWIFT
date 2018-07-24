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

def interpol_lambda(temp_list,cooling_rate,temp):
	#print(temp,temp_list[0],temp_list[len(temp_list)-1])
	if temp < temp_list[0]: 
		return[cooling_rate[0]]
	if temp > temp_list[len(temp_list)-1]: 
		return[cooling_rate[len(cooling_rate)-1]]
	j = 0
	while temp_list[j+1] < temp:
		j += 1
	cooling = cooling_rate[j]*(temp_list[j+1]-temp)/(temp_list[j+1]-temp_list[j]) + cooling_rate[j+1]*(temp-temp_list[j])/(temp_list[j+1]-temp_list[j])
	return cooling

def convert_u_to_temp_sol(u,rho):
	k_b = 1.38E-16 #boltzmann
	m_p = 1.67e-24 #proton mass
	gamma = 5.0/3.0
	n = 2.0*0.752*rho/m_p + 0.248*rho/(4.0*m_p)
	pressure = u*rho*(gamma - 1.0)
	temp = pressure/(k_b*n)
	#temp = pressure/(k_b*rho/(0.59*m_p))
	return temp

# File containing the total energy
stats_filename = "./energy.txt"

# First snapshot
snap_filename = "coolingBox_0000.hdf5"
snap_filename_mat = "../../../swiftsim_matthieu_cooling_v2/examples/CoolingBox/coolingBox_0000.hdf5"

# Some constants in cgs units
k_b = 1.38E-16 #boltzmann
m_p = 1.67e-24 #proton mass

# Initial conditions set in makeIC.py
T_init = 1.0e7

# Read the initial state of the gas
f = h5.File(snap_filename,'r')
f_mat = h5.File(snap_filename_mat,'r')
rho = np.mean(f["/PartType0/Density"])
rho_mat = np.mean(f["/PartType0/Density"])
pressure = np.mean(f["/PartType0/Pressure"])

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]
unit_vel = unit_length/unit_time
#hubble_param = 0.71
hubble_param = 1.0

unit_length = unit_length/hubble_param
unit_mass = unit_mass/hubble_param

yHe = 0.28
temp_0 = 1.0e7

rho = rho*unit_mass/(unit_length**3)
rho_mat = rho_mat*unit_mass/(unit_length**3)

# Read snapshots
nsnap = 501
npart = 4096
u_snapshots_cgs = zeros(nsnap)
u_part_snapshots_cgs = zeros((nsnap,npart))
t_snapshots_cgs = zeros(nsnap)
for i in range(nsnap):
    snap = h5.File("coolingBox_%0.4d.hdf5"%i,'r')
    u_part_snapshots_cgs[i][:] = snap["/PartType0/InternalEnergy"][:]*unit_length**2/(unit_time**2)
    t_snapshots_cgs[i] = snap["/Header"].attrs["Time"] * unit_time

# Read snapshots
#nsnap_mat = 25
#npart = 4096
#u_snapshots_cgs_mat = zeros(nsnap_mat)
#u_part_snapshots_cgs_mat = zeros((nsnap_mat,npart))
#t_snapshots_cgs_mat = zeros(nsnap_mat)
#for i in range(nsnap_mat):
#    snap = h5.File("../../../swiftsim_matthieu_cooling_v2/examples/CoolingBox/coolingBox_%0.4d.hdf5"%i,'r')
#    u_part_snapshots_cgs_mat[i][:] = snap["/PartType0/InternalEnergy"][:]*unit_length**2/(unit_time**2)
#    t_snapshots_cgs_mat[i] = snap["/Header"].attrs["Time"] * unit_time

# calculate analytic solution. Since cooling rate is constant,
# temperature decrease in linear in time.
# read Lambda and temperatures from table
temperature = []
cooling_rate = []
file_in = open('cooling_output.dat', 'r')
for line in file_in:
        data = line.split()
        temperature.append(float(data[0]))
        cooling_rate.append(-float(data[1]))

tfinal = t_snapshots_cgs[nsnap-1]
nt = 1e5
dt = tfinal/nt

t_sol = np.zeros(int(nt))
temp_sol = np.zeros(int(nt))
u_sol = np.zeros(int(nt))
lambda_sol = np.zeros(int(nt))
u = np.mean(u_part_snapshots_cgs[0,:])
temp_sol[0] = convert_u_to_temp_sol(u,rho)
#print(u,temp_sol[0])
u_sol[0] = u
for j in range(int(nt-1)):
	t_sol[j+1] = t_sol[j] + dt
	Lambda_net = interpol_lambda(temperature,cooling_rate,u_sol[j])
	nH = 0.702*rho/(m_p)
	#nH = 1.0e-4
	u_next = u - Lambda_net*nH**2/rho*dt
	temp_sol[j+1] = convert_u_to_temp_sol(u_next,rho)
	u_sol[j+1] = u_next
	#lambda_sol[j] = Lambda_net
	u = u_next
	
print(u, Lambda_net)

mean_temp = np.zeros(nsnap)
mean_u = np.zeros(nsnap)
for j in range(nsnap):
	mean_temp[j] = convert_u_to_temp_sol(np.mean(u_part_snapshots_cgs[j,:]),rho)
	mean_u[j] = np.mean(u_part_snapshots_cgs[j,:])
#mean_temp_mat = np.zeros(nsnap_mat)
#for j in range(nsnap_mat):
#	mean_temp_mat[j] = convert_u_to_temp_sol(np.mean(u_part_snapshots_cgs_mat[j,:]),rho_mat)

p = plt.figure(figsize = (6,6))
p1, = plt.loglog(t_sol,u_sol,linewidth = 0.5,color = 'k', marker = '*', markersize = 1.5,label = 'explicit ODE')
p2, = plt.loglog(t_snapshots_cgs,mean_u,linewidth = 0.5,color = 'r', marker = '*', markersize = 1.5,label = 'swift mean')
l = legend(handles = [p1,p2])
xlabel("${\\rm{Time~[s]}}$", labelpad=0, fontsize = 14)
ylabel("Internal energy ${\\rm{[erg \cdot g^{-1}]}}$", fontsize = 14)
title('$n_h = 1 \\rm{cm}^{-3}, z = 0$, zero metallicity,\n relative error: ' + "{0:.4f}".format( (u_sol[int(nt)-1] - mean_u[nsnap-1])/(u_sol[int(nt)-1])), fontsize = 16)
#plt.ylim([8.0e11,3.0e12])

savefig("energy.png", dpi=200)

print('Final internal energy ode, swift, error ' + str(u_sol[int(nt)-1]) + ' ' + str(mean_u[nsnap-1])  + ' ' + str( (u_sol[int(nt)-1] - mean_u[nsnap-1])/(u_sol[int(nt)-1])))
#print(u_sol[int(nt)-1], mean_u[nsnap-1], (u_sol[int(nt)-1] - mean_u[nsnap-1])/(u_sol[int(nt)-1]))
