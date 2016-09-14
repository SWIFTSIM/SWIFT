import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import sys

stats_filename = "./energy.txt"
snap_filename = "coolingBox_000.hdf5"
#plot_dir = "./"

#some constants in cgs units
k_b = 1.38E-16 #boltzmann
m_p = 1.67e-24 #proton mass
#initial conditions set in makeIC.py
rho = 3.2e3
P = 4.5e6
#n_H_cgs = 0.0001
gamma = 5./3.
T_init = 1.0e5

#Read the units parameters from the snapshot
f = h5.File(snap_filename,'r')
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]
parameters = f["Parameters"]
cooling_lambda = float(parameters.attrs["LambdaCooling:lambda_cgs"])
min_T = float(parameters.attrs["LambdaCooling:minimum_temperature"])
mu = float(parameters.attrs["LambdaCooling:mean_molecular_weight"])
X_H = float(parameters.attrs["LambdaCooling:hydrogen_mass_abundance"])

#get number of particles
header = f["Header"]
n_particles = header.attrs["NumPart_ThisFile"][0]

#read energy and time arrays
array = np.genfromtxt(stats_filename,skip_header = 1)
time = array[:,0]
kin_plus_therm = array[:,2]
radiated = array[:,6]
total_mass = array[:,1]

#ignore first row where there are just zeros
time = time[1:]
kin_plus_therm = kin_plus_therm[1:]
radiated = radiated[1:]
total_mass = total_mass[1:]

total_energy = kin_plus_therm + radiated
initial_energy = total_energy[0]
#conversions to cgs
rho_cgs = rho * unit_mass / (unit_length)**3
time_cgs = time * unit_time
initial_energy_cgs = initial_energy/total_mass[0] * unit_length**2 / (unit_time)**2 
n_H_cgs = X_H * rho_cgs / m_p

#find the energy floor
u_floor_cgs = k_b * min_T / (mu * m_p * (gamma - 1.))

#find analytic solution
analytic_time_cgs = np.linspace(0,time_cgs[-1],1000)
du_dt_cgs = -cooling_lambda * n_H_cgs**2 / rho_cgs
u_analytic_cgs = du_dt_cgs*analytic_time_cgs + initial_energy_cgs
cooling_time_cgs = initial_energy_cgs/(-du_dt_cgs)

for i in range(u_analytic_cgs.size):
    if u_analytic_cgs[i]<u_floor_cgs:
        u_analytic_cgs[i] = u_floor_cgs

#rescale analytic solution
u_analytic = u_analytic_cgs/initial_energy_cgs

#put time in units of cooling_time
time=time_cgs/cooling_time_cgs
analytic_time = analytic_time_cgs/cooling_time_cgs

#rescale (numerical) energy by initial energy
radiated /= initial_energy
kin_plus_therm /= initial_energy
total_energy = kin_plus_therm + radiated
plt.plot(time,kin_plus_therm,'kd',label = "Kinetic + thermal energy")
plt.plot(time,radiated,'bo',label = "Radiated energy")
plt.plot(time,total_energy,'g',label = "Total energy")
plt.plot(analytic_time,u_analytic,'r',lw = 2.0,label = "Analytic Solution")
#plt.plot(analytic_time,1-u_analytic,'k',lw = 2.0)
#plt.plot((cooling_time,cooling_time),(0,1),'b',label = "Cooling time")
#plt.plot((time[1]-time_cgs[0],time_cgs[1]-time_cgs[0]),(0,1),'m',label = "First output")
#plt.title(r"$n_H = %1.1e \, \mathrm{cm}^{-3}$" %n_H_cgs)
plt.xlabel("Time / cooling time")
plt.ylabel("Energy / Initial energy")
#plt.ylim(0,1.1)
plt.ylim(0.999,1.001)
#plt.xlim(0,min(10,time[-1]))
plt.legend(loc = "upper right")    
if (int(sys.argv[1])==0):
    plt.show()
else:
    plt.savefig(full_plot_filename,format = "png")
    plt.close()
