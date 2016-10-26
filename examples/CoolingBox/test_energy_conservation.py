import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import sys

stats_filename = "./energy.txt"
snap_filename = "coolingBox_000.hdf5"
#plot_dir = "./"
n_snaps = 41
time_end = 4.0
dt_snap = 0.1
#some constants in cgs units
k_b = 1.38E-16 #boltzmann
m_p = 1.67e-24 #proton mass
#initial conditions set in makeIC.py
rho = 4.8e3
P = 4.5e6
#n_H_cgs = 0.0001
gamma = 5./3.
T_init = 1.0e5

#find the sound speed

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
total_energy = array[:,2]
total_mass = array[:,1]

time = time[1:]
total_energy = total_energy[1:]
total_mass = total_mass[1:]

#conversions to cgs
rho_cgs = rho * unit_mass / (unit_length)**3
time_cgs = time * unit_time
u_init_cgs = total_energy[0]/(total_mass[0]) * unit_length**2 / (unit_time)**2 
n_H_cgs = X_H * rho_cgs / m_p

#find the sound speed in cgs
c_s = np.sqrt((gamma - 1.)*k_b*T_init/(mu*m_p))
#assume box size is unit length
sound_crossing_time = unit_length/c_s

print "Sound speed = %g cm/s" %c_s
print "Sound crossing time = %g s" %sound_crossing_time 
#find the energy floor
u_floor_cgs = k_b * min_T / (mu * m_p * (gamma - 1.))
#find analytic solution
analytic_time_cgs = np.linspace(time_cgs[0],time_cgs[-1],1000)
du_dt_cgs = -cooling_lambda * n_H_cgs**2 / rho_cgs
u_analytic = du_dt_cgs*(analytic_time_cgs - analytic_time_cgs[0]) + u_init_cgs
cooling_time = u_init_cgs/(-du_dt_cgs)

#put time in units of sound crossing time
time=time_cgs/sound_crossing_time
analytic_time = analytic_time_cgs/sound_crossing_time
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
plt.plot(time-time[0],total_energy,'k',label = "Numerical solution from energy.txt")
plt.plot(analytic_time-analytic_time[0],u_analytic,'r',lw = 2.0,label = "Analytic Solution")

#now get energies from the snapshots
snapshot_time = np.linspace(0,time_end,num = n_snaps)
snapshot_time = snapshot_time[1:]
snapshot_time_cgs = snapshot_time * unit_time
snapshot_time = snapshot_time_cgs/ sound_crossing_time
snapshot_time -= snapshot_time[0]
snapshot_energy = np.zeros(n_snaps)
for i in range(0,n_snaps):
    snap_filename = "coolingBox_%03d.hdf5" %i
    f = h5.File(snap_filename,'r')
    snapshot_internal_energy_array = np.array(f["PartType0/InternalEnergy"])
    total_internal_energy = np.sum(snapshot_internal_energy_array)
    velocity_array = np.array(f["PartType0/Velocities"])
    total_kinetic_energy = 0.5*np.sum(velocity_array**2)
    snapshot_energy[i] = total_internal_energy + total_kinetic_energy
snapshot_energy/=snapshot_energy[0]
snapshot_energy = snapshot_energy[1:]

plt.plot(snapshot_time,snapshot_energy,'bd',label = "Numerical solution from snapshots")

#plt.title(r"$n_H = %1.1e \, \mathrm{cm}^{-3}$" %n_H_cgs)
plt.xlabel("Time (sound crossing time)")
plt.ylabel("Energy/Initial energy")
plt.ylim(0.99,1.01)
#plt.xlim(0,min(10,time[-1]))
plt.legend(loc = "upper right")    
if (int(sys.argv[1])==0):
    plt.show()
else:
    plt.savefig(full_plot_filename,format = "png")
    plt.close()
