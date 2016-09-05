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
n_H_cgs = 0.0001
gamma = 5./3.
T_init = 1.0e5

#Read the units parameters from the snapshot
f = h5.File(snap_filename,'r')
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]
parameters = f["Parameters"]
cooling_lambda = float(parameters.attrs["LambdaCooling:lambda"])
min_T = float(parameters.attrs["LambdaCooling:minimum_temperature"])
mu = float(parameters.attrs["LambdaCooling:mean_molecular_weight"])
X_H = float(parameters.attrs["LambdaCooling:hydrogen_mass_abundance"])

#get number of particles
header = f["Header"]
n_particles = header.attrs["NumPart_ThisFile"][0]
#read energy and time arrays
array = np.genfromtxt(stats_filename,skip_header = 1)
time = array[:,0]
total_energy = array[:,2]
total_mass = array[:,1]

time = time[1:]
total_energy = total_energy[1:]
total_mass = total_mass[1:]

#conversions to cgs
rho_cgs = rho * unit_mass / (unit_length)**3
time_cgs = time * unit_time
u_init_cgs = total_energy[0]/(total_mass[0]) * unit_length**2 / (unit_time)**2 

#find the energy floor
print min_T
u_floor_cgs = k_b * min_T / (mu * m_p * (gamma - 1.))
#find analytic solution
analytic_time = np.linspace(time_cgs[0],time_cgs[-1],1000)
print time_cgs[1]
print analytic_time[1]
du_dt_cgs = -cooling_lambda * n_H_cgs**2 / rho_cgs
u_analytic = du_dt_cgs*(analytic_time - analytic_time[0]) + u_init_cgs
cooling_time = u_init_cgs/(-du_dt_cgs)
#rescale energy to initial energy
total_energy /= total_energy[0]
u_analytic /= u_init_cgs
u_floor_cgs /= u_init_cgs
# plot_title = r"$\Lambda \, = \, %1.1g \mathrm{erg}\mathrm{cm^3}\mathrm{s^{-1}} \, \, T_{init} = %1.1g\mathrm{K} \, \, T_{floor} = %1.1g\mathrm{K} \, \, n_H = %1.1g\mathrm{cm^{-3}}$" %(cooling_lambda,T_init,T_floor,n_H)
# plot_filename = "energy_plot_creasey_no_cooling_T_init_1p0e5_n_H_0p1.png"
#analytic_solution = np.zeros(n_snaps-1)
for i in range(u_analytic.size):
    if u_analytic[i]<u_floor_cgs:
        u_analytic[i] = u_floor_cgs
plt.plot(time_cgs,total_energy,'k',label = "Numerical solution")
plt.plot(analytic_time,u_analytic,'--r',lw = 2.0,label = "Analytic Solution")
plt.plot((cooling_time,cooling_time),(0,1),'b',label = "Cooling time")
plt.plot((time_cgs[0],time_cgs[0]),(0,1),'m',label = "First output")
plt.title(r"$n_H = %1.1e \, \mathrm{cm}^{-3}$" %n_H_cgs)
plt.xlabel("Time (seconds)")
plt.ylabel("Energy/Initial energy")
plt.ylim(0.999,1.001)
plt.xlim(0,min(10.0*cooling_time,time_cgs[-1]))
plt.legend(loc = "upper right")    
if (int(sys.argv[1])==0):
    plt.show()
else:
    plt.savefig(full_plot_filename,format = "png")
    plt.close()
